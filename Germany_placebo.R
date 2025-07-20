###### Placebo Test #####

rm(list=ls())
library(foreign)
library(Synth)
library(xtable)
library(kernlab)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(dplyr)

placebo_time = "1975-01-01"

# Load Data 
Germany <- read.csv("Germany.csv",header=TRUE)
Germany <- Germany %>% rename(Date = X)
Germany$Date <- as.Date(paste0(Germany$Date, "-01-01"))
colnames(Germany)[2:ncol(Germany)] <- c("Germany","Australia", "Austria", "Belgium", "Denmark","France","Greece","Italy","Japan","Netherlands","New Zealand","Norway","Portugal","Spain","Switzerland","UK","US")


h = 4
p = 2
Fh = h

T0 <- which(Germany$Date == placebo_time)

predata <- Germany[1:T0,]



# ========================= detrend all series ========================= 
## detrend function with OLS
lsq <- function(z){
  tt = length(z)
  y = z[(h+p):tt]
  X = embed(z[1:(tt-h)],p)
  y <- ifelse(is.na(y), 0, y)
  X <- ifelse(is.na(X), 0, X)
  df = data.frame(cbind(y,X))
  model = lm(y~.,data=df)
  return(model)
}

# function: Extract residuals and form a matrix
extract_residuals <- function(lm_list) {
  residuals_list <- lapply(lm_list, function(item) {
    return(residuals(item))
  })
  
  # Combine residuals into a matrix (each column corresponds to a model)
  residuals_matrix <- do.call(cbind, residuals_list)
  return(residuals_matrix)
}

trend_predict <- function(z, bhat) {
  tt = length(z)
  predictions <- numeric(h)
  for (i in 1:h){
    lags_h <- z[(tt-h+i):(tt-h-p+i+1)]
    x <- c(1,lags_h)
    predictions[i] <- sum(x*bhat)
  }
  return(predictions)
}

#-------- cyclical component of control units
list_detrended_Y0 <- apply(Germany[1:(T0+Fh),3:ncol(Germany)], 2, function(col) {
  result <- lsq(col)
  return(result)
})

detrended_donor <- extract_residuals(list_detrended_Y0)
detrended_Germany_pre <- lsq(predata[,2])$residuals

bhat_Germany = lsq(predata[,2])$coefficients
trend_Germany_pre <- tail(predata[,2],length(detrended_Germany_pre)) - detrended_Germany_pre
trend_Germany_post <- trend_predict(predata[,2], bhat_Germany)

trend_donor_pre <- tail(predata[,-c(1,2)],length(detrended_Germany_pre)) - detrended_donor[1:length(detrended_Germany_pre)]


trend_all_pre <- data.frame(Date = Germany$Date[(h+p):T0])
trend_all_pre$Germany <- trend_Germany_pre
trend_all_pre[,3:ncol(Germany)] <- trend_donor_pre
colnames(trend_all_pre)[2:ncol(Germany)] <- c("Germany","Australia", "Austria", "Belgium", "Denmark","France","Greece","Italy","Japan","Netherlands","New Zealand","Norway","Portugal","Spain","Switzerland","UK","US")


detrended_all_pre <- data.frame(Date = Germany$Date[(h+p):T0])
detrended_all_pre$Germany <- detrended_Germany_pre
detrended_all_pre[,3:ncol(Germany)] <- detrended_donor[1:length(detrended_Germany_pre),]
colnames(detrended_all_pre)[2:ncol(Germany)] <- c("Germany","Australia", "Austria", "Belgium", "Denmark","France","Greece","Italy","Japan","Netherlands","New Zealand","Norway","Portugal","Spain","Switzerland","UK","US")


# =================================================================#
# ============== sc/sbc restricted w >= 0 ==========================#
# =================================================================#

X1 <- as.matrix(predata$Germany)
X0 <- as.matrix(subset(predata[,2:ncol(predata)], select = -Germany))

V <- rep(1, T0)
synth.sc <- synth(data.prep.obj = NULL,
                  X1 = X1, X0 = X0,
                  Z0 = X0, Z1 = X1,
                  custom.v = V,
                  optimxmethod = c("Nelder-Mead", "BFGS"), 
                  genoud = FALSE, quadopt = "ipop", 
                  Margin.ipop = 0.0001, 
                  Sigf.ipop = 7, 
                  Bound.ipop = 6, 
                  verbose = FALSE)

weights.sc = synth.sc$solution.w
# SC estimates of Y1
donor_predict_df <- Germany[(h+p):(T0+Fh),3:ncol(Germany)]
Y1hat.SC <- tail(as.matrix(donor_predict_df),(T0+Fh-h-p+1))%*%weights.sc 
Germany_compare <- Germany$Germany[(h+p):(T0+Fh)]


X1 <- as.matrix(detrended_all_pre$Germany)
X0 <- as.matrix(subset(detrended_all_pre[,2:ncol(predata)], select = -Germany))

V <- rep(1, (T0-h-p+1))
synth.sbc <- synth(data.prep.obj = NULL,
                   X1 = X1, X0 = X0,
                   Z0 = X0, Z1 = X1,
                   custom.v = V,
                   optimxmethod = c("Nelder-Mead", "BFGS"), 
                   genoud = FALSE, quadopt = "ipop", 
                   Margin.ipop = 0.0001, 
                   Sigf.ipop = 7, 
                   Bound.ipop = 6, 
                   verbose = FALSE)

weights.sbc = synth.sbc$solution.w
cychat.sbc <- as.matrix(detrended_donor)%*%weights.sbc
trend_Germany <- c(trend_Germany_pre,trend_Germany_post)
Y1hat.sbc <- trend_Germany + cychat.sbc

#-------------------- Plots of SC/SBC restricted
time <- Germany$Date[(h+p):(T0+Fh)]
plot_data <- data.frame(
  time = rep(time, 3),
  value = c(Germany_compare, Y1hat.sbc, Y1hat.SC),
  group = rep(c("Germany", "SBC Germany", "SC Germany"), 
              each = length(time))
)

ggplot(plot_data, aes(x = time, y = value, color = group, linetype = group, size = group)) +
  geom_line() +
  geom_vline(xintercept = as.Date(placebo_time),
             linetype = "dotted", color = "black", size = 0.6) +
  scale_color_manual(values = c("Germany" = "black", 
                                "SBC Germany" = "black",
                                "SC Germany" = "grey")) +
  scale_linetype_manual(values = c("Germany" = "solid",
                                   "SBC Germany" = "dashed",
                                   "SC Germany" = "solid")) +
  scale_size_manual(values = c("Germany" = 1,
                               "SBC Germany" = 1,
                               "SC Germany" = 1)) +
  scale_x_date(
    breaks = seq.Date(from = as.Date("1965-01-01"),
                      to = max(plot_data$time, na.rm = TRUE),
                      by = "5 years"),  
    minor_breaks = seq.Date(from = min(plot_data$time, na.rm = TRUE),
                            to = max(plot_data$time, na.rm = TRUE),
                            by = "1 year"), 
    date_labels = "%Y",  
    expand = c(0.02, 0)  
  ) +
  labs(x = "", y = "Per Capita GDP (PPP, 2002 USD)") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),  
    legend.key.size = unit(1.2, "lines"),    
    legend.key.width = unit(1.2, "cm"),        
    legend.text = element_text(size = 14),      
    legend.margin = margin(2, 2, 2, 2),        
    legend.spacing = unit(0.2, "cm")          
  ) +
  guides(color = guide_legend(title = NULL),   
         linetype = guide_legend(title = NULL),
         size = guide_legend(title = NULL))


Y1diff.sbc <- Y1hat.sbc - Germany_compare
Y1diff.SC <- Y1hat.SC - Germany_compare
base <- rep(0,length(Germany_compare))
plot_diff <- data.frame(
  time = rep(time, 2),
  value = c(Y1diff.sbc, Y1diff.SC),
  group = rep(c("SBC Germany", "SC Germany"), 
              each = length(time))
)


ggplot(plot_diff, aes(x = time, y = value, color = group, linetype = group, size = group)) +
  geom_line() +
  geom_vline(xintercept = as.Date(placebo_time),
             linetype = "dotted", color = "black", size = 0.6) +
  geom_hline(yintercept = 0, color = "black", size = 1.2) +
  scale_color_manual(values = c(
    "SBC Germany" = "black",
    "SC Germany" = "grey")) +
  scale_linetype_manual(values = c(
    "SBC Germany" = "dashed",
    "SC Germany" = "solid")) +
  scale_size_manual(values = c(
    "SBC Germany" = 1,
    "SC Germany" = 1)) +
  scale_x_date(
    breaks = seq.Date(from = as.Date("1965-01-01"),
                      to = max(plot_data$time, na.rm = TRUE),
                      by = "5 years"),  
    minor_breaks = seq.Date(from = min(plot_data$time, na.rm = TRUE),
                            to = max(plot_data$time, na.rm = TRUE),
                            by = "1 year"),  
    date_labels = "%Y",  
    expand = c(0.02, 0) 
  ) +
  scale_y_continuous(limits = c(-2000, 4000)) +
  labs(x = "", y = "Per Capita GDP (PPP, 2002 USD)") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),  
    legend.key.size = unit(1.2, "lines"),      
    legend.key.width = unit(1.2, "cm"),        
    legend.text = element_text(size = 14),     
    legend.margin = margin(2, 2, 2, 2),        
    legend.spacing = unit(0.2, "cm")           
  ) +
  guides(color = guide_legend(title = NULL),   
         linetype = guide_legend(title = NULL),
         size = guide_legend(title = NULL))

