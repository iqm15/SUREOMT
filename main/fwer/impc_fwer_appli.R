library(dplyr)
library(OnlineSuperUnif)
library(onlineFDR)
library(DiscreteFDR)
library(lubridate)
library(rjson)
library(ggplot2)
library(tidyverse)
library(latex2exp)

# ---------------------------------------------------------------------------------------------------------
# impc_dataset <- read.csv("../impc/ReproducibleCode/Figure 4/SDofGenotypeEffect_processedData_categorical.csv")
# impc_dataset
# df <- select(impc_dataset, Male.Control.Count, Male_AbRateWT, Male.Mutant.Count, Male_AbRateKO, Female.Control.Count, Fem_AbRateWT, Female.Mutant.Count, Fem_AbRateKO, stage1__BH)
# df 
# write.csv(df, "my_impcdata.csv", row.names = TRUE)

# ------------------ load parameters from json file ------------------------------------------------------
param.list <- fromJSON(file = "../../config/fwer/impc_fwer_appli.json") 

# --------------------------------------------------------------------------------------------------------
impc_dataset <- read.csv(param.list$datasetpath)

# --------------------------------------------------------------------------------------------------------
# this function allows to recover the counts from the rates 
get_counts <- function(totals, rates) {
  if (length(totals) != length(rates)) { stop("The number of row must be the same")}
  Ab_counts <- numeric(length(totals)) # for abnormality counts
  N_counts <- numeric(length(totals)) # for normality counts
  for (i in 1:length(totals)) {
    Ab_counts[i] = rates[i] * totals[i]
    N_counts[i] = totals[i] - Ab_counts[i]
  }
  output <- list(Ab_counts = Ab_counts, N_counts = N_counts)
  return(output)
}


# -----------------------------------------------------------------------------------------------------------
Male_control <- get_counts(pull(impc_dataset, Male.Control.Count), pull(impc_dataset, Male_AbRateWT))
Male_mutant <- get_counts(pull(impc_dataset, Male.Mutant.Count), pull(impc_dataset, Male_AbRateKO))

Female_control <- get_counts(pull(impc_dataset, Female.Control.Count), pull(impc_dataset, Fem_AbRateWT))
Female_mutant <- get_counts(pull(impc_dataset, Female.Mutant.Count), pull(impc_dataset, Fem_AbRateKO))


# -----------------------------------------------------------------------------------------------------------
# creating dataset contiaining the contingency tables for each gene 
Male_df <- data.frame(cbind(Male_mutant$Ab_counts, Male_mutant$N_counts, 
                            Male_control$Ab_counts, Male_control$N_counts))
Fem_df <- data.frame(cbind(Female_mutant$Ab_counts, Female_mutant$N_counts, 
                           Female_control$Ab_counts, Female_control$N_counts))

# ------------------------------------------------------------------------------------------------------------
# write.csv(Male_df, "data/fwer/impc/impc_male_df.csv", row.names = TRUE)
# write.csv(Fem_df, "data/fwer/impc/impc_female_df.csv", row.names = TRUE)


# -------------------------------------------------------------------------------------------------------------
# compute the FET' p-values 
Male_test <- fisher.pvalues.support(Male_df, alternative = "greater", input = "noassoc")
Fem_test <- fisher.pvalues.support(Fem_df, alternative = "greater", input = "noassoc")


# ------------------------- PARAMETERS -------------------------------------------------------------------------
m_test = param.list$m_test # the number of genes to test must be between 1 and nrow(impc_dataset)
alpha = param.list$alpha
lambda = param.list$lambda
gamma <- gamma_sequence(param.list$gamma_type, m_test, param.list$q_1)
gamma_prime <- gamma_sequence(param.list$gamma_prime_type, m_test, param.list$q_2)


# ----------------------------- PROCEDURES ----------------------------------------------------------------------
# we perform the base OB, AOB and their rewarded counterparts rhoOB, rhoAOB for 
# female and male seqparetely 

OB_f <- OB(alpha, Fem_test$raw[1 : m_test], gamma)
OB_m <- OB(alpha, Male_test$raw[1 : m_test], gamma)

AOB_f <- AOB(alpha, Fem_test$raw[1 : m_test], gamma, lambda)
AOB_m <- AOB(alpha, Male_test$raw[1 : m_test], gamma, lambda)

rho_OB_f <- rho_OB(alpha, Fem_test$raw[1 : m_test], Fem_test$support[1 : m_test], gamma, gamma_prime)
rho_OB_m <- rho_OB(alpha, Male_test$raw[1 : m_test], Male_test$support[1 : m_test], gamma, gamma_prime)

rho_AOB_f <- rho_AOB(alpha, Fem_test$raw[1 : m_test], Fem_test$support[1 : m_test], gamma, lambda, gamma_prime)
rho_AOB_m <- rho_AOB(alpha, Male_test$raw[1 : m_test], Male_test$support[1 : m_test], gamma, lambda, gamma_prime)


# -------------------------- Create dataframes and save them ---------------------------------------------------------
# creates dataframes of critical values from the previous procedures runs

cvs_f <- data.frame(cbind(Fem_test$raw[1 : m_test], AOB_f$cv, OB_f$cv, rho_AOB_f$cv, rho_OB_f$cv))
names(cvs_f)[1] = "pvalues"
names(cvs_f)[2] = "AOB"
names(cvs_f)[3] = "OB"
names(cvs_f)[4] = "rho_AOB"
names(cvs_f)[5] = "rho_OB"
file_name = gsub(" " , "", paste("../../data/fwer/impc/", gsub(" ", "_", paste("impc_female_fwer_appli", now(), sep="_")), ".csv"))
write.csv(cvs_f, file_name, row.names = FALSE)

cvs_m <- data.frame(cbind(Male_test$raw[1 : m_test], AOB_m$cv, OB_m$cv, rho_AOB_m$cv, rho_OB_m$cv))
names(cvs_m)[1] = "pvalues"
names(cvs_m)[2] = "AOB"
names(cvs_m)[3] = "OB"
names(cvs_m)[4] = "rho_AOB"
names(cvs_m)[5] = "rho_OB"
file_name = gsub(" " , "", paste("../../data/fwer/impc/", gsub(" ", "_", paste("impc_male_fwer_appli", now(), sep="_")), ".csv"))
write.csv(cvs_m, file_name, row.names = FALSE)

# --------------------------- Create plots and save them -------------------------------------------------------------
m_plot = param.list$m_plot
# -log(-log) function used to transform the y-scale in the plots 
llog <- function(x) { 
  -log(-log(x))
}
# female 

# nb_rej_aob <- sum(cvs_f$pvalues <= cvs_f$AOB)
# nb_rej_ob <- sum(cvs_f$pvalues <= cvs_f$OB)
# nb_rej_rhoaob <- sum(cvs_f$pvalues <= cvs_f$rho_AOB)
# nb_rej_rhoob <- sum(cvs_f$pvalues <= cvs_f$rho_OB)
# nb_rej <- c(nb_rej_ob, nb_rej_rhoob, nb_rej_aob, nb_rej_rhoaob)
# nb_rej

data <- data.frame(seq(1, m_plot), llog(cvs_f$pvalues[1 : m_plot]), llog(cvs_f$OB[1 : m_plot]), llog(cvs_f$rho_OB[1 : m_plot]))
names(data)[1] = "index"
names(data)[2] = "pvalues"
names(data)[3] = "ob_cv"
names(data)[4] = "rhoob_cv"

data_rug_1 <- data.frame(seq(1, m_plot)[cvs_f$pvalues[1 : m_plot] <= cvs_f$rho_OB[1 : m_plot]])
names(data_rug_1)[1] <- "rhoobdisco"

data_rug_2 <- data.frame(seq(1, m_plot)[cvs_f$pvalues[1 : m_plot] <= cvs_f$OB[1 : m_plot]])
names(data_rug_2)[1] <- "obdisco"

plot_1_fem_loglog <- ggplot(data) +
  geom_point(aes_string("index", "pvalues"), size = 2, color = "black", alpha = 0.3) +
  geom_line(aes_string("index", "ob_cv"), lwd = 1, color = "darksalmon") +
  geom_line(aes_string("index", "rhoob_cv"), lwd = 1, color = "darkolivegreen") +
  geom_rug(data = data_rug_1, mapping = aes(x = rhoobdisco), color = "darkolivegreen", size = 1) +
  geom_rug(data = data_rug_2, mapping = aes(x = obdisco), color = "darksalmon", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "t", y = "p-values / critical values") +
  ylim(c(min(data$pvalues), 0)) +
  theme(axis.text = element_text(size = 12))

ggsave("../../figures/application/fwer/impc_app_fem_ob_rhoob.png", width = 16, height = 10, units = "cm")

data <- data.frame(seq(1, m_plot), llog(cvs_f$pvalues[1 : m_plot]), llog(cvs_f$AOB[1 : m_plot]), llog(cvs_f$rho_AOB[1 : m_plot]))
names(data)[1] = "index"
names(data)[2] = "pvalues"
names(data)[3] = "aob_cv"
names(data)[4] = "rhoaob_cv"

data_rug_1 <- data.frame(seq(1, m_plot)[cvs_f$pvalues[1 : m_plot] <= cvs_f$rho_AOB[1 : m_plot]])
names(data_rug_1)[1] <- "rhoaobdisco"

data_rug_2 <- data.frame(seq(1, m_plot)[cvs_f$pvalues[1 : m_plot] <= cvs_f$AOB[1 : m_plot]])
names(data_rug_2)[1] <- "aobdisco"

plot_2_fem_loglog <- ggplot(data) +
  geom_point(aes_string("index", "pvalues"), size = 2, color = "black", alpha = 0.3) +
  geom_line(aes_string("index", "aob_cv"), lwd = 1, color = "darksalmon") +
  geom_line(aes_string("index", "rhoaob_cv"), lwd = 1, color = "darkolivegreen") +
  geom_rug(data = data_rug_1, mapping = aes(x = rhoaobdisco), color = "darkolivegreen", size = 1) +
  geom_rug(data = data_rug_2, mapping = aes(x = aobdisco), color = "darksalmon", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "t", y = "p-values / critical values") +
  ylim(c(min(data$pvalues), 0)) +
  theme(axis.text = element_text(size = 12))

ggsave("../../figures/application/fwer/impc_app_fem_aob_rhoaob.png", width = 16, height = 10, units = "cm")

# male 

# nb_rej_aob <- sum(cvs_m$pvalues <= cvs_m$AOB)
# nb_rej_ob <- sum(cvs_m$pvalues <= cvs_m$OB)
# nb_rej_rhoaob <- sum(cvs_m$pvalues <= cvs_m$rho_AOB)
# nb_rej_rhoob <- sum(cvs_m$pvalues <= cvs_m$rho_OB)
# nb_rej <- c(nb_rej_ob, nb_rej_rhoob, nb_rej_aob, nb_rej_rhoaob)
# nb_rej

data <- data.frame(seq(1, m_plot), llog(cvs_m$pvalues[1 : m_plot]), llog(cvs_m$OB[1 : m_plot]), llog(cvs_m$rho_OB[1 : m_plot]))
names(data)[1] = "index"
names(data)[2] = "pvalues"
names(data)[3] = "ob_cv"
names(data)[4] = "rhoob_cv"

data_rug_1 <- data.frame(seq(1, m_plot)[cvs_m$pvalues[1 : m_plot] <= cvs_m$rho_OB[1 : m_plot]])
names(data_rug_1)[1] <- "rhoobdisco"

data_rug_2 <- data.frame(seq(1, m_plot)[cvs_m$pvalues[1 : m_plot] <= cvs_m$OB[1 : m_plot]])
names(data_rug_2)[1] <- "obdisco"

plot_1_male_loglog <- ggplot(data) +
  geom_point(aes_string("index", "pvalues"), size = 2, color = "black", alpha = 0.3) +
  geom_line(aes_string("index", "ob_cv"), lwd = 1, color = "darksalmon") +
  geom_line(aes_string("index", "rhoob_cv"), lwd = 1, color = "darkolivegreen") +
  geom_rug(data = data_rug_1, mapping = aes(x = rhoobdisco), color = "darkolivegreen", size = 1) +
  geom_rug(data = data_rug_2, mapping = aes(x = obdisco), color = "darksalmon", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "t", y = "p-values / critical values") +
  ylim(c(min(data$pvalues), 0)) +
  theme(axis.text = element_text(size = 12))

ggsave("../../figures/application/fwer/impc_app_male_ob_rhoob.png", width = 16, height = 10, units = "cm")

data <- data.frame(seq(1, m_plot), llog(cvs_m$pvalues[1 : m_plot]), llog(cvs_m$AOB[1 : m_plot]), llog(cvs_m$rho_AOB[1 : m_plot]))
names(data)[1] = "index"
names(data)[2] = "pvalues"
names(data)[3] = "aob_cv"
names(data)[4] = "rhoaob_cv"

data_rug_1 <- data.frame(seq(1, m_plot)[cvs_m$pvalues[1 : m_plot] <= cvs_m$rho_AOB[1 : m_plot]])
names(data_rug_1)[1] <- "rhoaobdisco"

data_rug_2 <- data.frame(seq(1, m_plot)[cvs_m$pvalues[1 : m_plot] <= cvs_m$AOB[1 : m_plot]])
names(data_rug_2)[1] <- "aobdisco"

plot_2_male_loglog <- ggplot(data) +
  geom_point(aes_string("index", "pvalues"), size = 2, color = "black", alpha = 0.3) +
  geom_line(aes_string("index", "aob_cv"), lwd = 1, color = "darksalmon") +
  geom_line(aes_string("index", "rhoaob_cv"), lwd = 1, color = "darkolivegreen") +
  geom_rug(data = data_rug_1, mapping = aes(x = rhoaobdisco), color = "darkolivegreen", size = 1) +
  geom_rug(data = data_rug_2, mapping = aes(x = aobdisco), color = "darksalmon", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "t", y = "p-values / critical values") +
  ylim(c(min(data$pvalues), 0)) +
  theme(axis.text = element_text(size = 12))

ggsave("../../figures/application/fwer/impc_app_male_aob_rhoaob.png", width = 16, height = 10, units = "cm")


# wealth plot (for male data here but could be done for female data also)
wealth_effective <- function(error_budget, cvs, stepf){
  # non-adaptive version
  m <- length(cvs)
  wealth.decrease <- sapply(1:m, FUN = function(i){stepf[[i]](cvs[i])})
  return(error_budget - cumsum(wealth.decrease))
}

step_m <- lapply(Male_test$support[1 : m_plot], function(x) stepfun(x, c(0, x)))

wealth.nominal_ob_m <- alpha - cumsum(cvs_m$OB)
wealth.effective_ob_m <- wealth_effective(alpha, cvs_m$OB, stepf_m)
wealth.nominal_rhoob_m <- alpha- cumsum(cvs_m$rho_OB)
wealth.effective_rhoob_m <- wealth_effective(alpha, cvs_m$rho_OB, stepf_m)

data_wealth_m <- data.frame(seq(1, length(wealth.nominal_ob_m)), wealth.nominal_ob_m, 
                            wealth.effective_ob_m, wealth.nominal_rhoob_m, wealth.effective_rhoob_m)
names(data_wealth_m)[1] = "index"
names(data_wealth_m)[2] = "nominal_wealth_ob"
names(data_wealth_m)[3] = "effective_wealth_ob"
names(data_wealth_m)[4] = "nominal_wealth_rhoob"
names(data_wealth_m)[5] = "effective_wealth_rhoob"

plot_male_wealth <- ggplot(data_wealth_m) +
  geom_line(aes_string("index", "nominal_wealth_ob"), lwd = 1, linetype = "dashed",  color = "darksalmon") +
  # geom_line(aes_string("index", "nominal_wealth_rhoob"), lwd = 1, linetype = "dashed", color = "darkolivegreen") +
  geom_line(aes_string("index", "effective_wealth_ob"), lwd = 1, color = "darksalmon") +
  geom_line(aes_string("index", "effective_wealth_rhoob"), lwd = 1, color = "darkolivegreen") +
  
  geom_hline(yintercept = llog(0.5), linetype = "dashed", color = "green", size = 0.5) +
  
  labs(tag = TeX("$\\lambda$")) +
  theme(plot.tag.position = c(0.03, 0.94)) +
  coord_cartesian( ylim = c(min(data$pvalues), llog(0.55)) , clip = "off") +
  labs(x = "t", y = "Wealth") +
  xlim(c(0, m_plot)) +
  # ylim(c(min(data_wealth_m$nominal_wealth_rhoob), alpha)) +
  theme(axis.text = element_text(size = 12))

ggsave("../../figures/application/fwer/male_wealth_ob_rhoob.png", width = 16, height = 10, units = "cm")
