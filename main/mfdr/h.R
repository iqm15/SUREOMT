setwd("/users/home/meah/phd_projects/SURE_OMT")
library(OnlineSuperUnif)
library(lubridate)
library(svMisc)
library(rjson)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(latex2exp)

# ------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config/mfdr/h.json")  
# -------------------------------------------------------------------------------

# nb of procedure to comapare
nb_proc = 2

# parameter of interest
non_nulls_prop <- seq(param.list$non_nulls_proportion$begin,
                      param.list$non_nulls_proportion$end,
                      by = param.list$non_nulls_proportion$by)  

kernel_bandwidths <- param.list$kernel_bandwidth  

# for collecting data on varying non_nulls_prop
means <- matrix(data = NA, nrow = nb_proc, ncol = length(non_nulls_prop), byrow = FALSE)  
sds <- matrix(data = NA, nrow = nb_proc, ncol = length(non_nulls_prop), byrow = FALSE)  
data_mfdr <- matrix(data = NA, nrow = nb_proc, ncol = length(non_nulls_prop), byrow = FALSE)  

# for collecting data on varying the kernel_bandwidth
means_ <- matrix(data=NA, nrow = nb_proc  * length(non_nulls_prop), ncol = length(kernel_bandwidths))
sds_ <- matrix(data=NA, nrow = nb_proc  * length(non_nulls_prop), ncol = length(kernel_bandwidths))
data_mfdr_ <- matrix(data=NA, nrow = nb_proc * length(non_nulls_prop), ncol = length(kernel_bandwidths))

# ----------------------------------------------------------------------------------------------------

for (j in 1:length(kernel_bandwidths)){
  for (i in 1:length(non_nulls_prop)) {
    data <- kernel_bandwidth_mfdr(param.list$alpha,
                                  param.list$w0,
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
    data_mfdr[, i] <- data$data_mfdr$mFDR_estimate
    
    means[, i] <- aggregate(as.numeric(data_power$ratio_true_discoveries),
                            by = list(data_power$procedure), FUN = mean)$x
    sds[, i] <- aggregate(as.numeric(data_power$ratio_true_discoveries),
                          by = list(data_power$procedure), FUN = sd)$x
    
  }
  means_[, j] = as.vector(means)
  sds_[, j] = as.vector(sds)
  data_mfdr_[, j] = as.vector(data_mfdr)
}
# ------------------------------------ create dataframe and save it  --------------------------------------------
procedure = rep(c("rhoALORD", "rhoLORD"), # alphabetical order 
                length.out = nb_proc * length(non_nulls_prop))  
param_interest.list = rep(non_nulls_prop, each = nb_proc)  

df_power <- data.frame(cbind(param_interest.list, procedure, means_, sds_))
df_mfdr <- data.frame(cbind(param_interest.list, procedure, data_mfdr_))

names(df_power)[1] = names(df_mfdr)[1] = param.list$param_interest
names(df_power)[2] = names(df_mfdr)[2] = "procedure"

for (i in 1:length(kernel_bandwidths)) {
  names(df_power)[i + 2] = gsub(" ", "_", paste("mean_bandwidth", kernel_bandwidths[i]))
  names(df_power)[i + 2 + length(kernel_bandwidths)] = gsub(" ", "_", paste("sd_bandwidth", kernel_bandwidths[i]))
  names(df_mfdr)[i + 2] = gsub(" ", "_", paste("mfdr_bandwidth", kernel_bandwidths[i]))
}

file_name_data_power = gsub(" " , "", paste("data/mfdr/h/", gsub(" ", "_", paste("data_power_bandwidth_vs_", param.list$param_interest, "study", now(), sep="_")), ".csv"))
file_name_data_mfdr = gsub(" " , "", paste("data/mfdr/h/", gsub(" ", "_", paste("data_mfdr_bandwidth_vs_", param.list$param_interest, "study", now(), sep="_")), ".csv"))
write.csv(df_power, file_name_data_power, row.names = FALSE)
write.csv(df_mfdr, file_name_data_mfdr, row.names = FALSE)

# --------------------------------- create plot and save it ----------------------------------------------------
plot_pow <- ggplot(df_power, aes_string(x = param.list$param_interest group = "procedure", linetype = "procedure")) +
  geom_line(aes(y = mean_bandwidth_1)) +
  geom_point(aes(y = mean_bandwidth_1, shape = procedure, color = "1"), size = 3) +
  
  # geom_line(aes(y = mean_bandwidth_2)) +
  # geom_point(aes(y = mean_bandwidth_2, shape = procedure, color = "2"), size = 3) +
  
  geom_line(aes(y = mean_bandwidth_5)) +
  geom_point(aes(y = mean_bandwidth_5, shape = procedure, color = "5"), size = 3) +
  
  geom_line(aes(y = mean_bandwidth_10)) +
  geom_point(aes(y = mean_bandwidth_10, shape = procedure, color = "10"), size = 3) +
  
  geom_line(aes(y = mean_bandwidth_50)) +
  geom_point(aes(y = mean_bandwidth_50, shape = procedure, color = "50"), size = 3) +
  
  geom_line(aes(y = mean_bandwidth_100)) +
  geom_point(aes(y = mean_bandwidth_100, shape = procedure, color = "100"), size = 3) +
  
  scale_color_discrete(breaks=c("1", "5", "10", "50", "100"), name = "Kernel bandwidth") +
  labs (x = TeX("$\\pi_A$") , y = "Power") +
  theme(text = element_text(size = 17))

saving_loc = "figures/simulation/mfdr/"
plot_name = gsub(" ", "", paste(saving_loc, param.list$param_interest, "kernel_bandwidth", ".png"))
ggsave(plot_name, plot = plot_pow)