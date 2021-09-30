library(OnlineSuperUnif)
library(lubridate)
library(svMisc)
library(rjson)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(latex2exp)

# ------------------ load parameters from json file -------------------------------------
param.list <- fromJSON(file = "../../config/fwer/h.json")  
# ---------------------------------------------------------------------------------------

# nb of procedure to compare
nb_proc = 2

# parameter of interest
non_nulls_prop <- seq(param.list$non_nulls_proportion$begin,
                      param.list$non_nulls_proportion$end,
                      by = param.list$non_nulls_proportion$by)  

kernel_bandwidths <- param.list$kernel_bandwidth  

# to collect data on varying non_nulls_prop
means <- matrix(data = NA, nrow = nb_proc, ncol = length(non_nulls_prop), byrow = FALSE)  
sds <- matrix(data = NA, nrow = nb_proc, ncol = length(non_nulls_prop), byrow = FALSE)  
data_fwer <- matrix(data = NA, nrow = nb_proc, ncol = length(non_nulls_prop), byrow = FALSE)  

# to collect data on varying the kernel_bandwidth
means_ <- matrix(data=NA, nrow = nb_proc  * length(non_nulls_prop), ncol = length(kernel_bandwidths))
sds_ <- matrix(data=NA, nrow = nb_proc  * length(non_nulls_prop), ncol = length(kernel_bandwidths))
data_fwer_ <- matrix(data=NA, nrow = nb_proc * length(non_nulls_prop), ncol = length(kernel_bandwidths))

# -----------------------------------------------------------------------------------------------------
for (j in 1:length(kernel_bandwidths)){
  for (i in 1:length(non_nulls_prop)) {
    data <- kernel_bandwidth_fwer(param.list$alpha,
                                  param.list$nb_run,
                                  param.list$N,
                                  param.list$m,
                                  non_nulls_prop[i],  
                                  param.list$p3,
                                  param.list$cluster_option,
                                  param.list$lambda,
                                  param.list$gamma_type,
                                  param.list$q_1,
                                  kernel_bandwidth = kernel_bandwidths[j])  
    
    data_power <- data$data_power
    data_fwer[, i] <- data$data_fwer$FWER_estimate
    
    means[, i] <- aggregate(as.numeric(data_power$ratio_true_discoveries),
                            by = list(data_power$procedure), FUN = mean)$x
    sds[, i] <- aggregate(as.numeric(data_power$ratio_true_discoveries),
                          by = list(data_power$procedure), FUN = sd)$x
    
  }
  means_[, j] = as.vector(means)
  sds_[, j] = as.vector(sds)
  data_fwer_[, j] = as.vector(data_fwer)
}

# -------------------------------- create data frames ans save it -------------------------------------
procedure = rep(c("rhoAOB", "rhoOB"),  # alphabetical order
                length.out = nb_proc * length(non_nulls_prop))  
param_interest.list = rep(non_nulls_prop, each = nb_proc)  

df_power <- data.frame(cbind(param_interest.list, procedure, means_, sds_))
df_fwer <- data.frame(cbind(param_interest.list, procedure, data_fwer_))

names(df_power)[1] = names(df_fwer)[1] = param.list$param_interest
names(df_power)[2] = names(df_fwer)[2] = "procedure"

for (i in 1:length(kernel_bandwidths)) {
  names(df_power)[i + 2] = gsub(" ", "_", paste("mean_bandwidth", kernel_bandwidths[i]))
  names(df_power)[i + 2 + length(kernel_bandwidths)] = gsub(" ", "_", paste("sd_bandwidth", kernel_bandwidths[i]))
  names(df_fwer)[i + 2] = gsub(" ", "_", paste("fwer_bandwidth", kernel_bandwidths[i]))
}

file_name_data_power = gsub(" " , "", paste( "../../data/fwer/h/", gsub(" ", "_", paste("data_power_fwer_bandwidth_vs_", param.list$param_interest, "study", now(), sep="_")), ".csv"))
file_name_data_fwer = gsub(" " , "", paste("../../data/fwer/h/", gsub(" ", "_", paste("data_fwer_bandwidth_vs_", param.list$param_interest, "study", now(), sep="_")), ".csv"))
write.csv(df_power, file_name_data_power, row.names = FALSE)
write.csv(df_fwer, file_name_data_fwer, row.names = FALSE)

# ------------------------------- create plot and save it ----------------------------------------------

data_power <- read.csv(file_name_data_power)

plot_pow <- ggplot(data_power, aes_string(x = param.list$param_interest,
                                       group = "procedure", 
                                       linetype = "procedure")) +
  geom_line(aes(y = mean_bandwidth_10)) +
  geom_point(aes(y = mean_bandwidth_10, shape = procedure, color = "10"), size = 3) +
  
  geom_line(aes(y = mean_bandwidth_30)) +
  geom_point(aes(y = mean_bandwidth_30, shape = procedure, color = "30"), size = 3) +
  
  geom_line(aes(y = mean_bandwidth_50)) +
  geom_point(aes(y = mean_bandwidth_50, shape = procedure, color = "50"), size = 3) +
  
  geom_line(aes(y = mean_bandwidth_70)) +
  geom_point(aes(y = mean_bandwidth_70, shape = procedure, color = "70"), size = 3) +
  
  geom_line(aes(y = mean_bandwidth_90)) +
  geom_point(aes(y = mean_bandwidth_90, shape = procedure, color = "90"), size = 3) +
  
  geom_line(aes(y = mean_bandwidth_100)) +
  geom_point(aes(y = mean_bandwidth_100, shape = procedure, color = "100"), size = 3) +
  
  scale_color_discrete(breaks=c("10", "30", "50", "70", "90", "100"), name = "Kernel bandwidth") +
  labs (x = TeX("$\\pi_A$") , y = "Power") +
  theme(text = element_text(size = 17))

saving_loc = "../../figures/simulation/fwer/"
plot_name = gsub(" ", "", paste(saving_loc, param.list$param_interest, ".png"))
ggsave(plot_name, plot = plot_pow)
