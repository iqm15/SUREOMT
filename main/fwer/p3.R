setwd("/users/home/meah/phd_projects/online_superunif_mt/code/sure_omt_test")
library(OnlineSuperUnif)
library(lubridate)
library(svMisc)
library(rjson)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(latex2exp)

#------------------ load parameters from yaml file -----------------------------
param.list <- fromJSON(file = "config/fwer/p3.json") ##
#-------------------------------------------------------------------------------
# number of procedures to compare
nb_proc = 5

# parameter of interest
p3 <- seq(param.list$p3$begin, param.list$p3$end, by = param.list$p3$by) ##

# to collect data to plot the power
means <- matrix(data = NA, nrow = nb_proc, ncol = length(p3), byrow = FALSE) ##
sds <- matrix(data = NA, nrow = nb_proc, ncol = length(p3), byrow = FALSE) ##

# to collect data to plot the fwer
data_fwer <- matrix(data = NA, nrow = nb_proc, ncol = length(p3), byrow = FALSE) ##

for (i in 1:length(p3)) {   ##
  data <- compare_fwer_proc(param.list$alpha, param.list$nb_run,
                            param.list$N, param.list$m,
                            param.list$non_nulls_prop, p3[i],   ##
                            param.list$cluster_option, param.list$lambda,
                            param.list$gamma_type, param.list$q_1,
                            param.list$gamma_prime_type, param.list$q_2)

  data_power <- data$data_power
  data_fwer[, i] <- data$data_fwer$fwer_estimate

  means[, i] <- aggregate(as.numeric(data_power$ratio_true_discoveries),
                          by = list(data_power$procedure), FUN = mean)$x
  sds[, i] <- aggregate(as.numeric(data_power$ratio_true_discoveries),
                        by = list(data_power$procedure), FUN = sd)$x

}
means = as.vector(means)
sds = as.vector(sds)
data_fwer = as.vector(data_fwer)


# -------------------- create data frame and save it ----------------------------

procedure = rep(c("ADDIS-spending", "AOB", "OB", "rhoAOB", "rhoOB"), # alphabetical order
                length.out = nb_proc * length(p3))   
param_interest.list = rep(p3, each = nb_proc)   

df_power <- data.frame(cbind(param_interest.list, procedure, means, sds))
df_fwer <- data.frame(cbind(param_interest.list, procedure, data_fwer))

names(df_power)[1] = names(df_fwer)[1] = param.list$param_interest
names(df_power)[2] = names(df_fwer)[2] = "procedure"

names(df_power)[3] = "mean"
names(df_power)[4] = "sd"
names(df_fwer)[3] = "fwer"

file_name_data_power = gsub(" " , "", paste("data/fwer/p3/", gsub(" ", "_", paste("data_power", param.list$param_interest, "study", now(), sep="_")), ".csv"))
file_name_data_fwer = gsub(" " , "", paste("data/fwer/p3/", gsub(" ", "_", paste("data_fwer", param.list$param_interest, "study", now(), sep="_")), ".csv"))
write.csv(df_power, file_name_data_power, row.names = FALSE)
write.csv(df_fwer, file_name_data_fwer, row.names = FALSE)

# -------------------- make plot and save it -----------------------------------------
# Power plot
plot_pow <- ggplot(df_power, 
                   aes_string(x = param.list$param_interest, 
                              y = "mean", 
                              group = "procedure", 
                              color = "procedure", 
                              linetype = "procedure", 
                              shape = "procedure") 
) +
  scale_color_brewer(palette = "Dark2") +
  geom_line(size = 1) +
  geom_point(size = 4) +
  scale_shape_manual(values = seq(0, 15)) + 
  labs (x = TeX("$p_3$") , y = "Power") +
  theme(text = element_text(size = 17))

# FWER plot
plot_fwer <- ggplot(df_fwer, 
                    aes_string(x = param.list$param_interest, 
                               y = "fwer", 
                               group = "procedure", 
                               color = "procedure", 
                               linetype = "procedure", 
                               shape = "procedure")) +
  scale_color_brewer(palette = "Dark2") +
  geom_line(size = 1) +
  geom_point(size = 4) +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "red", size = 0.5) +
  scale_shape_manual(values = seq(0, 15)) +
  labs (x = TeX("$p_3$"), y = "FWER") +
  theme(text = element_text(size = 17))

# arrange both plots in one figure and save it 
figure <- ggarrange(plot_pow, plot_fwer, ncol = 2, nrow = 1, common.legend = TRUE)
saving_loc = "figures/simulation/fwer/"
plot_name = gsub(" ", "", paste(saving_loc, param.list$param_interest, ".png"))
ggsave(plot_name, plot = figure)
