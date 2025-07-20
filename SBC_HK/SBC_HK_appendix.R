# install.packages("readxl")  # 安装包（如果未安装）

## SC, SBC with w = 1, 9 donors

library(readxl)
library(dplyr)

library(ggplot2)
library(tidyr)
library(quadprog)

library(foreign)
library(Synth)
library(xtable)
library(kernlab)
library(gridExtra)

# load in real GDP data, 1960 - 2010

# Get a list of all CSV files in the folder
csv_files <- list.files(path = folder_path, pattern = "\\.csv$", full.names = TRUE)

# Load all CSV files into a list
csv_list <- lapply(csv_files, read.csv)

# Optionally, assign names to the list elements based on file names
names(csv_list) <- basename(csv_files)

# View the first dataset
head(csv_list[[1]])

# Ensure all datasets have Date column formatted as Date type
csv_list <- lapply(csv_list, function(df) {
  df[[1]] <- as.Date(df[[1]])  # Assuming first column is Date
  df
})   

treat_time = "1997-01-01"

date_col <- csv_list[[1]][, 1]  

gdp_df <- data.frame(Date = date_col)

for (i in seq_along(csv_list)) {
  country_name <- tools::file_path_sans_ext(names(csv_list)[i])  # 去除.csv后缀
  gdp_df[[country_name]] <- csv_list[[i]][, 2]  # 添加第二列数据
}

str(gdp_df)
head(gdp_df)
T = length(gdp_df$Date)
T0 <- which(gdp_df$Date == treat_time)

h = 4
p = 2
Fh = h


# ========================= Plot HK and 24 donor countries ========================= 

donor_df <- gdp_df 
ind_HK = which(colnames(donor_df) == "HK")
cols <- names(donor_df)
donor_df <- donor_df[, c(1, ind_HK, setdiff(seq_along(cols), c(1, ind_HK)))]

donor_df <- donor_df %>%
  select(Date, HK, China, Indonesia,
         Japan, Korea, Malaysia, Philippines, Singapore, Thailand, US)


# Convert data from wide to long format for ggplot
df_long <- donor_df %>%
  pivot_longer(cols = -Date, names_to = "Country", values_to = "Value")

df_long <- df_long %>%
  mutate(Color = ifelse(Country == "HK", "HK", "Donor Countries"))



regular_breaks <- seq(
  from = min(as.Date("1960-01-01")),
  to = max(df_long$Date),
  by = "10 years"
)


ggplot(df_long, aes(x = Date, y = Value, group = Country, color = Color)) +
  geom_line(data = df_long[df_long$Country != "HK", ], size = 0.5) +
  geom_line(data = df_long[df_long$Country == "HK", ], size = 1) +
  geom_vline(xintercept = as.Date(treat_time), linetype = "dotted", color = "black") +
  scale_x_date(
    breaks = regular_breaks,
    labels = function(x) format(x, "%Y")
  ) +
  scale_color_manual(values = c("HK" = "black", "Donor Countries" = "grey")) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "PPP, 2010 USD", color = "Legend") +
  theme(legend.position = "none")

# ======================================SC with non-negative weights ======================================
predata <- donor_df[1:T0,]
X1 <- as.matrix(predata$HK)
X0 <- as.matrix(subset(predata[,2:ncol(predata)], select = -HK))

V <- rep(1, T0)
synth.sc <- synth(data.prep.obj = NULL,
                  X1 = X1, X0 = X0,
                  Z0 = X0, Z1 = X1,
                  custom.v = V,
                  optimxmethod = c("Nelder-Mead", "BFGS"), 
                  genoud = FALSE, quadopt = "ipop", 
                  Margin.ipop = 0.0005, 
                  Sigf.ipop = 7, 
                  Bound.ipop = 6, 
                  verbose = FALSE)

weights.sc = synth.sc$solution.w

# SC estimates 
donor_predict_df <- donor_df[(h+p):(T0+Fh),3:ncol(donor_df)]
Y1hat.SC <- as.matrix(donor_predict_df)%*%weights.sc
HK_compare <- donor_df$HK[(h+p):(T0+Fh)]

# ====================================== SBC restricted weights ======================================
## detrend function with OLS
lsq <- function(z){
  tt = length(z)
  y = z[(h+p):tt]
  X = embed(z[1:(tt-h)],p)
  #X = cbind(1,X)
  y <- ifelse(is.na(y), 0, y)
  X <- ifelse(is.na(X), 0, X)
  df = data.frame(cbind(y,X))
  model = lm(y~.,data=df)
  return(model)
}

# function: Extract residuals and form a matrix
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
  predictions_temp <- c(z,numeric(Fh))
  
  
  for (i in (tt+1):(tt+Fh)){
    lags <- predictions_temp[(i-h):(i-h-p+1)]
    # lags_h <- z[(tt-h+i):(tt-h-p+i+1)]
    x <- c(1,lags)
    predictions_temp[i] <- sum(x*bhat)
  }
  predictions = predictions_temp[(tt+1):(tt+Fh)]
  
  return(predictions)
}

time<- donor_df$Date[(h+p):(T0+Fh)]
detrend_df <- donor_df[1:(T0+Fh),2:ncol(donor_df)]
donor_detrend_df <-  subset(detrend_df, select = -HK) # donor detrend 1:T0+Fh
HK_detrend_df <- detrend_df$HK[1:T0] # HK detrend 1:T0


list_detrended_control <- apply(donor_detrend_df, 2, function(col) {
  result <- lsq(col)
  return(result)
})
donor_detrended <- extract_residuals(list_detrended_control) # h+p : T0+Fh
HK_detrended <- lsq(HK_detrend_df)$residuals 
donor_detrended_pre <- donor_detrended[1:length(HK_detrended),]

HK_trend_pre = tail(HK_detrend_df,length(HK_detrended)) - HK_detrended
HK_trend_coef = lsq(HK_detrend_df)$coefficients

# predict HK trend after treatment with linear projection
HK_trend_post <- trend_predict(HK_detrend_df, HK_trend_coef)
HK_trend_hat <- c(HK_trend_pre,HK_trend_post)
# predict UK cyclical after treatment with synthetic
synth.SBC <- synth(data.prep.obj = NULL,
                   X1 = as.matrix(HK_detrended), X0 = as.matrix(donor_detrended_pre),
                   Z0 = as.matrix(donor_detrended_pre), Z1 = as.matrix(HK_detrended),
                   custom.v = rep(1,length(HK_detrended)),
                   optimxmethod = c("Nelder-Mead", "BFGS"), 
                   genoud = FALSE, quadopt = "ipop", 
                   Margin.ipop = 5e-04, 
                   Sigf.ipop = 7, 
                   Bound.ipop = 6, 
                   verbose = FALSE)
weights.SBC <- synth.SBC$solution.w
HK_cyc_hat <- donor_detrended%*%weights.SBC
Y1hat.SBC <- HK_trend_hat+HK_cyc_hat

detrended_all_pre <- donor_df[(h + p):T0, 1, drop = FALSE]
detrended_all_pre <-cbind(detrended_all_pre, HK_detrended, donor_detrended_pre)
names(detrended_all_pre)[2] <- "HK"

trend_all_pre <- detrend_df[(h+p):T0,] - detrended_all_pre[,2:ncol(detrended_all_pre)]
trend_all_pre <- cbind(detrended_all_pre[,1], trend_all_pre)
names(trend_all_pre)[2] <- "HK"
names(trend_all_pre)[1] <- "Date"

df_long <- trend_all_pre %>%
  pivot_longer(cols = -Date, names_to = "Country", values_to = "Value")

df_long <- df_long %>%
  mutate(Color = ifelse(Country == "HK", "HK", "Donor Countries"))

regular_breaks <- seq(
  from = min(as.Date("1965-01-01")),
  to = max(df_long$Date),
  by = "5 years"
)

ggplot(df_long, aes(x = Date, y = Value, group = Country, color = Color)) +
  geom_line(data = df_long[df_long$Country != "HK", ], size = 0.5) +
  geom_line(data = df_long[df_long$Country == "HK", ], size = 1) +
  geom_vline(xintercept = as.Date(treat_time), linetype = "dotted", color = "black") +
  scale_x_date(
    breaks = regular_breaks,
    labels = function(x) format(x, "%Y")
  ) +
  scale_color_manual(values = c("HK" = "black", "Donor Countries" = "grey")) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "PPP, 2010 USD", color = "Legend") +
  theme(legend.position = "none")

#-------------------- Plots of SC/SBC with non-negative weights
time <- gdp_df$Date[(h+p):(T0+Fh)]
plot_data <- data.frame(
  time = rep(time, 3),
  value = c(HK_compare, Y1hat.SBC, Y1hat.SC),
  group = rep(c("HK", "SBC HK", "SC HK"), 
              each = length(time))
)

ggplot(plot_data, aes(x = time, y = value, color = group, linetype = group, size = group)) +
  geom_line() +
  geom_vline(xintercept = as.Date(treat_time),
             linetype = "dotted", color = "black", size = 0.6) +
  scale_color_manual(values = c("HK" = "black", 
                                "SBC HK" = "black",
                                "SC HK" = "grey")) +
  scale_linetype_manual(values = c("HK" = "solid",
                                   "SBC HK" = "dashed",
                                   "SC HK" = "solid")) +
  scale_size_manual(values = c("HK" = 1,
                               "SBC HK" = 1,
                               "SC HK" = 1)) +
  scale_x_date(
    breaks = seq.Date(from = as.Date("1967-01-01"),
                      to = max(plot_data$time, na.rm = TRUE),
                      by = "5 years"),  
    minor_breaks = seq.Date(from = min(plot_data$time, na.rm = TRUE),
                            to = max(plot_data$time, na.rm = TRUE),
                            by = "1 year"),  
    date_labels = "%Y",  
    expand = c(0.02, 0)  
  ) +
  labs(x = "", y = "Per Capita GDP (PPP, 2010 USD)") +
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


# ========================= Plot all cyclical before T0 =========================


df_long <- detrended_all_pre %>%
  pivot_longer(cols = -Date, names_to = "Country", values_to = "Value")

df_long <- df_long %>%
  mutate(Color = ifelse(Country == "HK", "HK", "Donor Countries"))

regular_breaks <- seq(
  from = min(as.Date("1965-01-01")),
  to = max(df_long$Date),
  by = "5 years"
)

ggplot(df_long, aes(x = Date, y = Value, color = Color)) +
  geom_line(data = df_long[df_long$Country != "HK", ],
            aes(x = Date, y = Value, group = Country), size = 0.5) +
  geom_line(data = df_long[df_long$Country == "HK", ],
            aes(x = Date, y = Value, group = Country), size = 1) +
  geom_vline(xintercept = as.Date(treat_time), linetype = "dotted", color = "black") +
  scale_x_date(
    breaks = regular_breaks,
    labels = function(x) format(x, "%Y")
  ) +
  scale_color_manual(values = c("HK" = "black", "Donor Countries" = "grey")) +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "PPP, 2010 USD", color = "Legend") +
  theme(legend.position = "none")

#-------------------- Plots of SC/SBC with non-negative weights
plot_data <- data.frame(
  time = rep(time, 3),
  value = c(HK_compare, Y1hat.SBC, Y1hat.SC),
  group = rep(c("HK", "SBC HK", "SC HK"),
              each = length(time))
)

ggplot(plot_data, aes(x = time, y = value, color = group, linetype = group, size = group)) +
  geom_line() +
  geom_vline(xintercept = as.Date(treat_time),
             linetype = "dotted", color = "black", size = 0.6) +
  scale_color_manual(values = c("HK" = "black",
                                "SBC HK" = "black",
                                "SC HK" = "grey")) +
  scale_linetype_manual(values = c("HK" = "solid",
                                   "SBC HK" = "dashed",
                                   "SC HK" = "solid")) +
  scale_size_manual(values = c("HK" = 1,
                               "SBC HK" = 1,
                               "SC HK" = 1)) +
  scale_x_date(
    breaks = seq.Date(from = as.Date("1960-01-01"),
                      to = max(plot_data$time, na.rm = TRUE),
                      by = "5 years"),  
    minor_breaks = seq.Date(from = min(plot_data$time, na.rm = TRUE),
                            to = max(plot_data$time, na.rm = TRUE),
                            by = "1 year"),  
    date_labels = "%Y",  
    expand = c(0.02, 0)  
  ) +
  labs(x = "", y = "Per Capita GDP (PPP, 2010 USD)") +
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

# ====================================== SC/SBC with signed weights sum to one ======================================
predata <- donor_df[1:T0,]
X1 <- as.matrix(predata$HK)
X0 <- as.matrix(subset(predata[,2:ncol(predata)], select = -HK))

D <- t(X0) %*% X0  # X'X
d <- t(X0) %*% X1  # X'y

A <- matrix(1, nrow = 1, ncol = ncol(X0))  
b <- 1  

sol <- solve.QP(D, d, t(A), b, meq = 1)  
beta_SC_w1 <- sol$solution
Y1hat_SC.w1 <- as.matrix(donor_predict_df)%*%beta_SC_w1

# SBC 
X1 = as.matrix(HK_detrended)
X0 = as.matrix(donor_detrended_pre)
D <- t(X0) %*% X0  # X'X
d <- t(X0) %*% X1  # X'y

A <- matrix(1, nrow = 1, ncol = ncol(X0)) 
b <- 1 

sol <- solve.QP(D, d, t(A), b, meq = 1)  
beta_SBC_w1 <- sol$solution
HK_cyc.w1 <- as.matrix(donor_detrended)%*%beta_SBC_w1
Y1hat_SBC.w1 <- HK_trend_hat+HK_cyc.w1


#-------------------- Plots of SC/SBC signed weights sum to one
# Define time variable

plot_data <- data.frame(
  time = rep(time, 3),
  value = c(HK_compare, Y1hat_SBC.w1, Y1hat_SC.w1),
  group = rep(c("HK", "SBC HK", "SC HK"), 
              each = length(time))
)


ggplot(plot_data, aes(x = time, y = value, color = group, linetype = group, size = group)) +
  geom_line() +
  geom_vline(xintercept = as.Date(treat_time),
             linetype = "dotted", color = "black", size = 0.6) +
  scale_color_manual(values = c("HK" = "black", 
                                "SBC HK" = "black",
                                "SC HK" = "grey")) +
  scale_linetype_manual(values = c("HK" = "solid",
                                   "SBC HK" = "dashed",
                                   "SC HK" = "solid")) +
  scale_size_manual(values = c("HK" = 1,
                               "SBC HK" = 1,
                               "SC HK" = 1)) +
  scale_x_date(
    breaks = seq.Date(from = as.Date("1960-01-01"),
                      to = max(plot_data$time, na.rm = TRUE),
                      by = "5 years"),  
    minor_breaks = seq.Date(from = min(plot_data$time, na.rm = TRUE),
                            to = max(plot_data$time, na.rm = TRUE),
                            by = "1 year"),  
    date_labels = "%Y",  
    expand = c(0.02, 0) 
  ) +
  labs(x = "", y = "Per Capita GDP (PPP, 2010 USD)") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = c(0.05, 0.9),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),  # 图例键背景
    legend.key.size = unit(1.2, "lines"),    
    legend.key.width = unit(1.2, "cm"),      
    legend.text = element_text(size = 14),     
    legend.margin = margin(0.5, 0.5, 0.5, 0.5),        
    legend.spacing = unit(0.2, "cm")           
  ) +
  guides(color = guide_legend(title = NULL),   
         linetype = guide_legend(title = NULL),
         size = guide_legend(title = NULL))

growth.SBC <- c(NA, diff(Y1hat.SBC) / Y1hat.SBC[-nrow(Y1hat.SBC)])
postgrowth.SBC <- tail(growth.SBC,Fh)

growth <- c(NA, diff(log(Y1hat.SBC)))
postgrowth <- tail(growth,Fh)

growth.HK <- c(NA, diff(log(HK_compare)))
postgrowth.HK <- tail(growth.HK,Fh)

tail(Y1hat.SBC,Fh+5)

# ======================== placebo test  ======================== #
placebo_date = "1987-01-01"
T0 <- which(gdp_df$Date == placebo_date )
predata <- donor_df[1:T0,]
X1 <- as.matrix(predata$HK)
X0 <- as.matrix(subset(predata[,2:ncol(predata)], select = -HK))

V <- rep(1, T0)
synth.sc <- synth(data.prep.obj = NULL,
                  X1 = X1, X0 = X0,
                  Z0 = X0, Z1 = X1,
                  custom.v = V,
                  optimxmethod = c("Nelder-Mead", "BFGS"), 
                  genoud = FALSE, quadopt = "ipop", 
                  Margin.ipop = 0.0005, 
                  Sigf.ipop = 7, 
                  Bound.ipop = 6, 
                  verbose = FALSE)

weights.sc = synth.sc$solution.w

# SC estimates 
donor_predict_df <- donor_df[(h+p):(T0+Fh),3:ncol(donor_df)]
Y1hat.SC <- as.matrix(donor_predict_df)%*%weights.sc
HK_compare <- donor_df$HK[(h+p):(T0+Fh)]


time<- donor_df$Date[(h+p):(T0+Fh)]
detrend_df <- donor_df[1:(T0+Fh),2:ncol(donor_df)]
donor_detrend_df <-  subset(detrend_df, select = -HK) # donor detrend 1:T0+Fh
HK_detrend_df <- detrend_df$HK[1:T0] # HK detrend 1:T0


list_detrended_control <- apply(donor_detrend_df, 2, function(col) {
  result <- lsq(col)
  return(result)
})
donor_detrended <- extract_residuals(list_detrended_control) # h+p : T0+Fh
HK_detrended <- lsq(HK_detrend_df)$residuals 
donor_detrended_pre <- donor_detrended[1:length(HK_detrended),]

HK_trend_pre = tail(HK_detrend_df,length(HK_detrended)) - HK_detrended
HK_trend_coef = lsq(HK_detrend_df)$coefficients

# predict HK trend after treatment with linear projection
HK_trend_post <- trend_predict(HK_detrend_df, HK_trend_coef)
HK_trend_hat <- c(HK_trend_pre,HK_trend_post)
# predict UK cyclical after treatment with synthetic
synth.SBC <- synth(data.prep.obj = NULL,
                   X1 = as.matrix(HK_detrended), X0 = as.matrix(donor_detrended_pre),
                   Z0 = as.matrix(donor_detrended_pre), Z1 = as.matrix(HK_detrended),
                   custom.v = rep(1,length(HK_detrended)),
                   optimxmethod = c("Nelder-Mead", "BFGS"), 
                   genoud = FALSE, quadopt = "ipop", 
                   Margin.ipop = 5e-04, 
                   Sigf.ipop = 7, 
                   Bound.ipop = 6, 
                   verbose = FALSE)
weights.SBC <- synth.SBC$solution.w
HK_cyc_hat <- donor_detrended%*%weights.SBC
Y1hat.SBC <- HK_trend_hat+HK_cyc_hat

detrended_all_pre <- donor_df[(h + p):T0, 1, drop = FALSE]
detrended_all_pre <-cbind(detrended_all_pre, HK_detrended, donor_detrended_pre)
names(detrended_all_pre)[2] <- "HK"

trend_all_pre <- detrend_df[(h+p):T0,] - detrended_all_pre[,2:ncol(detrended_all_pre)]
trend_all_pre <- cbind(detrended_all_pre[,1], trend_all_pre)
names(trend_all_pre)[2] <- "HK"
names(trend_all_pre)[1] <- "Date"

df_long <- trend_all_pre %>%
  pivot_longer(cols = -Date, names_to = "Country", values_to = "Value")

df_long <- df_long %>%
  mutate(Color = ifelse(Country == "HK", "HK", "Donor Countries"))

regular_breaks <- seq(
  from = min(as.Date("1965-01-01")),
  to = max(df_long$Date),
  by = "5 years"
)

ggplot(df_long, aes(x = Date, y = Value, group = Country, color = Color)) +
  geom_line(data = df_long[df_long$Country != "HK", ], size = 0.5) +
  geom_line(data = df_long[df_long$Country == "HK", ], size = 1) +
  geom_vline(xintercept = as.Date(treat_time), linetype = "dotted", color = "black") +
  scale_x_date(
    breaks = regular_breaks,
    labels = function(x) format(x, "%Y")
  ) +
  scale_color_manual(values = c("HK" = "black", "Donor Countries" = "grey")) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_minimal(base_size = 16) +
  labs(x = "", y = "PPP, 2010 USD", color = "Legend") +
  theme(legend.position = "none")

#-------------------- Plots of SC/SBC with non-negative weights
time <- gdp_df$Date[(h+p):(T0+Fh)]

plot_data <- data.frame(
  time = rep(time, 3),
  value = c(HK_compare, Y1hat.SBC, Y1hat.SC),
  group = rep(c("HK", "SBC HK", "SC HK"),
              each = length(time))
)

ggplot(plot_data, aes(x = time, y = value, color = group, linetype = group, size = group)) +
  geom_line() +
  geom_vline(xintercept = as.Date(treat_time),
             linetype = "dotted", color = "black", size = 0.6) +
  scale_color_manual(values = c("HK" = "black",
                                "SBC HK" = "black",
                                "SC HK" = "grey")) +
  scale_linetype_manual(values = c("HK" = "solid",
                                   "SBC HK" = "dashed",
                                   "SC HK" = "solid")) +
  scale_size_manual(values = c("HK" = 1,
                               "SBC HK" = 1,
                               "SC HK" = 1)) +
  scale_x_date(
    breaks = seq.Date(from = as.Date("1960-01-01"),
                      to = max(plot_data$time, na.rm = TRUE),
                      by = "5 years"),  
    minor_breaks = seq.Date(from = min(plot_data$time, na.rm = TRUE),
                            to = max(plot_data$time, na.rm = TRUE),
                            by = "1 year"),  
    date_labels = "%Y",  
    expand = c(0.02, 0)  
  ) +
  labs(x = "", y = "Per Capita GDP (PPP, 2010 USD)") +
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

# ====================================== SC/SBC with signed weights sum to one ======================================
predata <- donor_df[1:T0,]
X1 <- as.matrix(predata$HK)
X0 <- as.matrix(subset(predata[,2:ncol(predata)], select = -HK))

D <- t(X0) %*% X0  # X'X
d <- t(X0) %*% X1  # X'y

A <- matrix(1, nrow = 1, ncol = ncol(X0))  
b <- 1  

sol <- solve.QP(D, d, t(A), b, meq = 1)  
beta_SC_w1 <- sol$solution
Y1hat_SC.w1 <- as.matrix(donor_predict_df)%*%beta_SC_w1

# SBC 
X1 = as.matrix(HK_detrended)
X0 = as.matrix(donor_detrended_pre)
D <- t(X0) %*% X0  # X'X
d <- t(X0) %*% X1  # X'y

A <- matrix(1, nrow = 1, ncol = ncol(X0)) 
b <- 1 

sol <- solve.QP(D, d, t(A), b, meq = 1)  
beta_SBC_w1 <- sol$solution
HK_cyc.w1 <- as.matrix(donor_detrended)%*%beta_SBC_w1
Y1hat_SBC.w1 <- HK_trend_hat+HK_cyc.w1


#-------------------- Plots of SC/SBC signed weights sum to one
# Define time variable

plot_data <- data.frame(
  time = rep(time, 3),
  value = c(HK_compare, Y1hat_SBC.w1, Y1hat_SC.w1),
  group = rep(c("HK", "SBC HK", "SC HK"), 
              each = length(time))
)


ggplot(plot_data, aes(x = time, y = value, color = group, linetype = group, size = group)) +
  geom_line() +
  geom_vline(xintercept = as.Date(treat_time),
             linetype = "dotted", color = "black", size = 0.6) +
  scale_color_manual(values = c("HK" = "black", 
                                "SBC HK" = "black",
                                "SC HK" = "grey")) +
  scale_linetype_manual(values = c("HK" = "solid",
                                   "SBC HK" = "dashed",
                                   "SC HK" = "solid")) +
  scale_size_manual(values = c("HK" = 1,
                               "SBC HK" = 1,
                               "SC HK" = 1)) +
  scale_x_date(
    breaks = seq.Date(from = as.Date("1960-01-01"),
                      to = max(plot_data$time, na.rm = TRUE),
                      by = "5 years"),  
    minor_breaks = seq.Date(from = min(plot_data$time, na.rm = TRUE),
                            to = max(plot_data$time, na.rm = TRUE),
                            by = "1 year"),  
    date_labels = "%Y",  
    expand = c(0.02, 0) 
  ) +
  labs(x = "", y = "Per Capita GDP (PPP, 2010 USD)") +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = c(0.05, 0.9),
    legend.justification = c(0, 1),
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),  # 图例键背景
    legend.key.size = unit(1.2, "lines"),    
    legend.key.width = unit(1.2, "cm"),      
    legend.text = element_text(size = 14),     
    legend.margin = margin(0.5, 0.5, 0.5, 0.5),        
    legend.spacing = unit(0.2, "cm")           
  ) +
  guides(color = guide_legend(title = NULL),   
         linetype = guide_legend(title = NULL),
         size = guide_legend(title = NULL))

growth.SBC <- c(NA, diff(Y1hat.SBC) / Y1hat.SBC[-nrow(Y1hat.SBC)])
postgrowth.SBC <- tail(growth.SBC,Fh)

growth <- c(NA, diff(log(Y1hat.SBC)))
postgrowth <- tail(growth,Fh)

growth.HK <- c(NA, diff(log(HK_compare)))
postgrowth.HK <- tail(growth.HK,Fh)

tail(Y1hat.SBC,Fh+5)
