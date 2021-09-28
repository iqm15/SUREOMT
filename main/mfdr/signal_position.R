setwd("/users/home/meah/phd_projects/online_superunif_mt/code/sure_omt_test")
library(OnlineSuperUnif)
library(lubridate)
library(svMisc)
library(rjson)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(latex2exp)

# ------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config/mfdr/signal_position.json")  
# -------------------------------------------------------------------------------
# number of procedures to compare
nb_proc = 6

# parameter of interest
cluster_option <- param.list$cluster_option  

# to collect data to plot the power
means <- matrix(data = NA, nrow = nb_proc, ncol = length(cluster_option), byrow = FALSE)  
sds <- matrix(data = NA, nrow = nb_proc, ncol = length(cluster_option), byrow = FALSE)  

# to collect data to plot the mfdr
data_mfdr <- matrix(data = NA, nrow = nb_proc, ncol = length(cluster_option), byrow = FALSE)  

# ------------------------------------------------------------------------------------------
for (i in 1:length(cluster_option)) {    
  data <- compare_mfdr_proc(param.list$alpha, param.list$w0,
                            param.list$nb_run,
                            param.list$N, param.list$m,
                            param.list$non_nulls_proportion, param.list$p3,
                            cluster_option[i], param.list$lambda,  
                            param.list$gamma_type, param.list$q_1,
                            param.list$gamma_prime_type, param.list$q_2)
  
  data_power <- data$data_power
  data_mfdr[, i] <- data$data_mfdr$mFDR_estimate
  
  means[, i] <- aggregate(as.numeric(data_power$ratio_true_discoveries),
                          by = list(data_power$procedure), FUN = mean)$x
  sds[, i] <- aggregate(as.numeric(data_power$ratio_true_discoveries),
                        by = list(data_power$procedure), FUN = sd)$x
  
}
means = as.vector(means)
sds = as.vector(sds)
data_mfdr = as.vector(data_mfdr)


# -------------------- create data frame and save it ----------------------------

procedure = rep(c("ADDIS", "ALORD", "LORD", "rhoALORD", "rhoLORD", "SAFFRON"), # alphabetical order 
                length.out = nb_proc * length(cluster_option))  
param_interest.list = rep(cluster_option, each = nb_proc)  

df_power <- data.frame(cbind(param_interest.list, procedure, means, sds))
df_mfdr <- data.frame(cbind(param_interest.list, procedure, data_mfdr))

names(df_power)[1] = names(df_mfdr)[1] = param.list$param_interest
names(df_power)[2] = names(df_mfdr)[2] = "procedure"

names(df_power)[3] = "mean"
names(df_power)[4] = "sd"
names(df_mfdr)[3] = "mfdr"

file_name_data_power = gsub(" " , "", paste("data/mfdr/signalpos/", gsub(" ", "_", paste("data_power", param.list$param_interest, "study", now(), sep="_")), ".csv"))
file_name_data_mfdr = gsub(" " , "", paste("data/mfdr/signalpos/", gsub(" ", "_", paste("data_mfdr", param.list$param_interest, "study", now(), sep="_")), ".csv"))
write.csv(df_power, file_name_data_power, row.names = FALSE)
write.csv(df_mfdr, file_name_data_mfdr, row.names = FALSE)

# -------------------- make plot and save it -------------------------------------
# to the x-ticks sorted in alphabetical order
nb_rows <- nrow(df_power)
for (i in 1:nb_rows) {
  if ( df_power[i, 1] == "begin") { df_power[i, 1] = df_mfdr[i, 1] = "a"}
  if ( df_power[i, 1] == "begin_middle") { df_power[i, 1] = df_mfdr[i, 1] = "b"}
  if ( df_power[i, 1] == "begin_end") { df_power[i, 1] = df_mfdr[i, 1] = "c"}
  if ( df_power[i, 1] == "middle_end") { df_power[i, 1] = df_mfdr[i, 1] = "d"}
  if ( df_power[i, 1] == "end") { df_power[i, 1] = df_mfdr[i, 1] = "e"}
  if ( df_power[i, 1] == "no_cluster_shuffle") { df_power[i, 1] = df_mfdr[i, 1] = "f"}
}
# Power plot 
plot_pow <- ggplot(df_power, 
                   aes_string(x = param.list$param_interest, 
                              y = "mean", 
                              group = "procedure", 
                              color = "procedure", 
                              linetype = "procedure", 
                              shape = "procedure") 
) +
  scale_color_brewer(palette = "Paired") +
  scale_x_discrete(labels=c("a" = "B",
                            "b" = "B.M" ,
                            "c" = "B.E",
                            "d" = "M.E",
                            "e" = "E",
                            "f" = "Random")) +
  
  geom_point(size = 4) +
  scale_shape_manual(values = seq(0, 15)) + 
  labs (x = "Signal position", y = "Power") +
  theme(text = element_text(size = 17))

# mFDR plot
plot_mfdr <- ggplot(df_mfdr, 
                    aes_string(x = param.list$param_interest, 
                               y = "mfdr", 
                               group = "procedure", 
                               color = "procedure", 
                               linetype = "procedure", 
                               shape = "procedure")) +
  scale_color_brewer(palette = "Paired") +
  scale_x_discrete(labels=c("a" = "B",
                            "b" = "B.M" ,
                            "c" = "B.E",
                            "d" = "M.E",
                            "e" = "E",
                            "f" = "Random")) +

  geom_point(size = 4) +
  geom_hline(yintercept = param.list$alpha, linetype = "dashed", color = "red", size = 0.5) +
  scale_shape_manual(values = seq(0, 15)) + 
  labs (x = "Signal Position", y = "mFDR") +
  theme(text = element_text(size = 17))

# arrange both plots in one figure and save it 
figure <- ggarrange(plot_pow, plot_mfdr, ncol = 2, nrow = 1, common.legend = TRUE)
saving_loc = "figures/simulation/mfdr/"
plot_name = gsub(" ", "", paste(saving_loc, param.list$param_interest, ".png"))
ggsave(plot_name, plot = figure)
