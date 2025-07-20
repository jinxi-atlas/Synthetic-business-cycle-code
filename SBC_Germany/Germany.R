rm(list=ls())
library(foreign)
library(Synth)
library(xtable)
library(kernlab)
library(ggplot2)
library(tidyr)
library(gridExtra)
library(dplyr)
library(ggpattern)
library(quadprog)

# Load Data 
Germany <- read.csv("Germany.csv",header=TRUE)
Germany <- Germany %>% rename(Date = X)
Germany$Date <- as.Date(paste0(Germany$Date, "-01-01"))

h = 4 # detrending horizon
p = 2 # number of lags used in detrending
Fh = h # forecast horizon after treatment

treat_time = "1990-01-01"
T0 <- which(Germany$Date == treat_time)

predata <- Germany[1:T0,] # pre-treatment data

# ========================= Plot Germany and all donor countries gdp ========================= 


# Convert data from wide to long format for ggplot
df_long <- Germany %>%
  pivot_longer(cols = -Date, names_to = "Country", values_to = "Value")

# Define colors: Germany in blue, others in grey
df_long <- df_long %>%
  mutate(Color = ifelse(Country == "Germany", "Germany", "Donor Countries"))

regular_breaks <- seq(
  from = min(as.Date("1960-01-01")),
  to = max(df_long$Date),
  by = "10 years"
)

ggplot(df_long, aes(x = Date, y = Value, group = Country, color = Color)) +
  geom_line(data = df_long[df_long$Country != "Germany", ], size = 0.5) +
  geom_line(data = df_long[df_long$Country == "Germany", ], size = 1) +
  geom_vline(xintercept = as.Date(treat_time), linetype = "dotted", color = "black") +
  scale_x_date(
    breaks = regular_breaks,
    labels = function(x) format(x, "%Y")
  ) +
  scale_color_manual(values = c("Germany" = "black", "Donor Countries" = "grey")) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "PPP, 2002 USD", color = "Legend") +
  theme(legend.position = "none")

# ========================= detrend all series ========================= 
# detrending with linear projection
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

# function: Extract residuals 
extract_residuals <- function(lm_list) {
  # Use lapply to extract residuals
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

#-------- calculate cyclical components
# list_detrended_Y0 <- apply(Germany[1:(T0+Fh),3:ncol(Germany)], 2, function(col) {
#   result <- lsq(col)
#   return(result)
# })

list_detrended_Y0 <- apply(Germany[,3:ncol(Germany)], 2, function(col) {
  result <- lsq(col)
  return(result)
})

detrended_donor <- extract_residuals(list_detrended_Y0) # cyclical components of donors
detrended_donor <- detrended_donor[1:(T0+Fh- h-p+1),]
detrended_Germany_pre <- lsq(predata[,2])$residuals # cyclical components of donors before treatment

bhat_Germany = lsq(predata[,2])$coefficients # detrending coefficients of Germany

# pre-treatment trend of Germany
trend_Germany_pre <- tail(predata[,2],length(detrended_Germany_pre)) - detrended_Germany_pre

# estimated post-treatment trend of Germany using pre-treatment values
trend_Germany_post <- trend_predict(predata[,2], bhat_Germany)

# pre-treatment trends of donors
trend_donor_pre <- tail(predata[,-c(1,2)],length(detrended_Germany_pre)) - detrended_donor[1:length(detrended_Germany_pre),]

# pre-treatment trends of all countries
trend_all_pre <- data.frame(Date = Germany$Date[(h+p):T0])

# pre-treatment trend in dataframe
trend_all_pre$Germany <- trend_Germany_pre
trend_all_pre[,3:ncol(Germany)] <- trend_donor_pre
colnames(trend_all_pre)[2:ncol(Germany)] <- c("Germany","Australia", "Austria", "Belgium", "Denmark","France","Greece","Italy","Japan","Netherlands","New Zealand","Norway","Portugal","Spain","Switzerland","UK","US")

# pre-treatment cyclical components in dataframe
detrended_all_pre <- data.frame(Date = Germany$Date[(h+p):T0])
detrended_all_pre$Germany <- detrended_Germany_pre
detrended_all_pre[,3:ncol(Germany)] <- detrended_donor[1:length(detrended_Germany_pre),]
colnames(detrended_all_pre)[2:ncol(Germany)] <- c("Germany","Australia", "Austria", "Belgium", "Denmark","France","Greece","Italy","Japan","Netherlands","New Zealand","Norway","Portugal","Spain","Switzerland","UK","US")



# ========================= Plot all trends before T0 ========================= 

df_long <- trend_all_pre %>%
  pivot_longer(cols = -Date, names_to = "Country", values_to = "Value")

# Define colors: Germany in blue, others in grey
df_long <- df_long %>%
  mutate(Color = ifelse(Country == "Germany", "Germany", "Donor Countries"))

regular_breaks <- seq(
  from = min(as.Date("1965-01-01")),
  to = max(df_long$Date),
  by = "5 years"
)

ggplot(df_long, aes(x = Date, y = Value, group = Country, color = Color)) +
  geom_line(data = df_long[df_long$Country != "Germany", ], size = 0.5) +
  geom_line(data = df_long[df_long$Country == "Germany", ], size = 1) +
  geom_vline(xintercept = as.Date(treat_time), linetype = "dotted", color = "black") +
  scale_x_date(
    breaks = regular_breaks,
    labels = function(x) format(x, "%Y")
  ) +
  scale_color_manual(values = c("Germany" = "black", "Donor Countries" = "grey")) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "PPP, 2002 USD", color = "Legend") +
  theme(legend.position = "none")


# ========================= Plot all cyclical before T0 ========================= 

df_long <- detrended_all_pre %>%
  pivot_longer(cols = -Date, names_to = "Country", values_to = "Value")

df_long <- df_long %>%
  mutate(Color = ifelse(Country == "Germany", "Germany", "Donor Countries"))

regular_breaks <- seq(
  from = min(as.Date("1965-01-01")),
  to = max(df_long$Date),
  by = "5 years"
)


ggplot(df_long, aes(x = Date, y = Value, color = Color)) +
  geom_line(data = df_long[df_long$Country != "Germany", ],
            aes(x = Date, y = Value, group = Country), size = 0.5) +
  geom_line(data = df_long[df_long$Country == "Germany", ],
            aes(x = Date, y = Value, group = Country), size = 1) +
  geom_vline(xintercept = as.Date(treat_time), linetype = "dotted", color = "black") +
  scale_x_date(
    breaks = regular_breaks,
    labels = function(x) format(x, "%Y")
  ) +
  scale_color_manual(values = c("Germany" = "black", "Donor Countries" = "grey")) +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "PPP, 2002 USD", color = "Legend") +
  theme(legend.position = "none")


# =================================================================#
# ============== sc/sbc restricted ,non-negative weights ==========#
# =================================================================#

# define treated - level of Germany
X1 <- as.matrix(predata$Germany) 

# define controls - level of donors
X0 <- as.matrix(subset(predata[,2:ncol(predata)], select = -Germany)) 

V <- rep(1, T0)

# perform synthetic control
synth.sc <- synth(data.prep.obj = NULL,
                  X1 = X1, X0 = X0,
                  Z0 = X0, Z1 = X1,
                  custom.v = V,
                  optimxmethod = c("Nelder-Mead", "BFGS"), 
                  genoud = FALSE, quadopt = "ipop", 
                  Margin.ipop = 0.0005, 
                  Sigf.ipop = 5, 
                  Bound.ipop = 10, 
                  verbose = FALSE)


# synthetic control weights
weights.sc = synth.sc$solution.w

# use donor to form synthetic control estimates
donor_predict_df <- Germany[(h+p):(T0+Fh),3:ncol(Germany)]
Y1hat.SC <- tail(as.matrix(donor_predict_df),(T0+Fh-h-p+1))%*%weights.sc 



# define treated - detrended Germany 
X1 <- as.matrix(detrended_all_pre$Germany)

# define control - detrended donors
X0 <- as.matrix(subset(detrended_all_pre[,2:ncol(predata)], select = -Germany))

# perform synthetic business cycle
V <- rep(1, (T0-h-p+1))
synth.sbc <- synth(data.prep.obj = NULL,
                   X1 = X1, X0 = X0,
                   Z0 = X0, Z1 = X1,
                   custom.v = V,
                   optimxmethod = c("Nelder-Mead", "BFGS"), 
                   genoud = FALSE, quadopt = "ipop", 
                   Margin.ipop = 5e-04, 
                   Sigf.ipop = 5, 
                   Bound.ipop = 10, 
                   verbose = FALSE)

# synthetic business cycle weights
weights.sbc = synth.sbc$solution.w

# form SBC estimates for Germany
cychat.sbc <- as.matrix(detrended_donor)%*%weights.sbc
trend_Germany <- c(trend_Germany_pre,trend_Germany_post)
Y1hat.sbc <- trend_Germany + cychat.sbc

#==================== Plots of SC/SBC non-negative weights - at levels  #
time <- Germany$Date[(h+p):(T0+Fh)]
Germany_compare <- Germany$Germany[(h+p):(T0+Fh)]

plot_data <- data.frame(
  time = rep(time, 3),
  value = c(Germany_compare, Y1hat.sbc, Y1hat.SC),
  group = rep(c("Germany", "SBC Germany", "SC Germany"), 
              each = length(time))
)

ggplot(plot_data, aes(x = time, y = value, color = group, linetype = group, size = group)) +
  geom_line() +
  geom_vline(xintercept = as.Date(treat_time),
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

#==================== Plots of SC/SBC non-negative weights - difference from actual Germany GDP  #

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
  geom_vline(xintercept = as.Date(treat_time),
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
    legend.position = c(0.05, 0.9),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(), 
    legend.key.size = unit(1.2, "lines"),      
    legend.key.width = unit(1.2, "cm"),        
    legend.text = element_text(size = 14),      
    legend.margin = margin(0.5, 0.5, 0.5, 0.5),       
    legend.spacing = unit(0.2, "cm")         
  ) +
  guides(color = guide_legend(title = NULL),
         linetype = guide_legend(title = NULL),
         size = guide_legend(title = NULL))

#==================== weights comparison of SC/SBC weights, w restricted ========================= #


synth.trendSC <- synth(data.prep.obj = NULL,
                       X1 = as.matrix(trend_Germany_pre), X0 = as.matrix(trend_donor_pre),
                       Z0 = as.matrix(trend_donor_pre), Z1 = as.matrix(trend_Germany_pre),
                       custom.v = rep(1,length(trend_Germany_pre)),
                       optimxmethod = c("Nelder-Mead", "BFGS"), 
                       genoud = FALSE, quadopt = "ipop", 
                       Margin.ipop = 5e-04, 
                       Sigf.ipop = 5, 
                       Bound.ipop = 10, 
                       verbose = FALSE)

weight.trendSC <- synth.trendSC$solution.w



# ---------------- plot weights ---------------------
barplot(
  rbind(c(weight.trendSC),c(weights.sc), c(weights.sbc)),
  beside = TRUE,
  col = c("white","black","grey"),
  names.arg = colnames(Germany[,3:ncol(Germany)]),
  legend.text = c("SC trend","SC", "SBC"),
  ylab = "weights",
  las = 2,         
  cex.names = 0.7,  
  args.legend = list(cex = 0.7)  
)


# ============== calculate errors ==========================#
error_sc <- tail(Y1diff.SC,Fh)
loss_sc <- sum(error_sc)/Fh


error_sbc <- tail(Y1diff.sbc,Fh)
loss_sbc <- sum(error_sbc)/Fh

