# -----------------------------------------------------------------------------
setwd("/users/home/meah/phd_projects/online_superunif_mt/code/sure_omt_test/")
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
param.list <- fromJSON(file = "config/mfdr/impc_mfdr_appli.json") 

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
w0 <- alpha * (1 - lambda)
gamma <- gamma_sequence(param.list$gamma_type, m_test, param.list$q_1)
gamma_prime <- gamma_sequence(param.list$gamma_prime_type, m_test, param.list$q_2) 

# ----------------------------- PROCEDURES ----------------------------------------------------------------------
# we perform the base LORD, ALORD and their rewarded counterparts rhoLORD, rhoALORD for 
# female and male seqparetely 
lord_f <- lord_OnlineSuperUnif(alpha, w0, Fem_test$raw[1 : m_test], gamma)
lord_m <- lord_OnlineSuperUnif(alpha, w0, Male_test$raw[1 : m_test], gamma)

# sanity check to verify that my implementation of LORD and the one 
# form the onlineFDR package provide the same results 
lord_onlinefdr_f <- LORD(Fem_test$raw[1 : m_test], alpha, gamma, w0 = w0)
lord_onlinefdr_m <- LORD(Male_test$raw[1 : m_test], alpha, gamma, w0 = w0)

sum((lord_onlinefdr_f$alphai - lord_f$cv)**2)
sum((lord_onlinefdr_m$alphai_m - lord_m$cv)**2)


# sanity check to verify that my implementation of SAFFRON and the one 
# form the onlineFDR package provide the same results
saffron_f <- saffron_OnlineSuperUnif(alpha, w0, Fem_test$raw[1 : m_test], gamma, lambda, capping = TRUE)
saffron_m <- saffron_OnlineSuperUnif(alpha, w0, Male_test$raw[1 : m_test], gamma, lambda, capping = TRUE)

saffron_onlinefdr_f <- SAFFRON(Fem_test$raw[1 : m_test], alpha, gamma, w0, lambda)
saffron_onlinefdr_m <- SAFFRON(Male_test$raw[1 : m_test], alpha, gamma, w0, lambda)

sum((saffron_f$cv - saffron_onlinefdr_f$alphai)**2)
sum((saffron_m$cv - saffron_onlinefdr_m$alphai)**2)


# for the comparison between procedures we use ALORD (and not SAFFRON)
alord_f <- saffron_OnlineSuperUnif(alpha, w0, Fem_test$raw[1 : m_test], gamma, lambda, capping = TRUE)
alord_m <- saffron_OnlineSuperUnif(alpha, w0, Male_test$raw[1 : m_test], gamma, lambda, capping = TRUE)

rholord_f <- rho_LORD(alpha, w0, Fem_test$raw[1 : m_test], Fem_test$support[1 : m_test], gamma, gamma_prime)
rholord_m <- rho_LORD(alpha, w0, Male_test$raw[1 : m_test], Male_test$support[1 : m_test], gamma, gamma_prime)

rhoalord_f <- rho_ALORD(alpha, w0, Fem_test$raw[1 : m_test], Fem_test$support[1 : m_test], gamma, lambda, gamma_prime)
rhoalord_m <- rho_ALORD(alpha, w0, Male_test$raw[1 : m_test], Male_test$support[1 : m_test], gamma, lambda, gamma_prime)



# -------------------------- Create dataframes and save them ---------------------------------------------------------
# creates dataframes of critical values from the previous procedures runs

cvs_f <- data.frame(cbind(Fem_test$raw[1 : m_test], alord_f$cv, lord_f$cv, rhoalord_f$cv, rholord_f$cv))
names(cvs_f)[1] = "pvalues"
names(cvs_f)[2] = "alord"
names(cvs_f)[3] = "lord"
names(cvs_f)[4] = "rho_alord"
names(cvs_f)[5] = "rho_lord"
file_name = gsub(" " , "", paste("data/mfdr/impc/", gsub(" ", "_", paste("impc_female_mfdr_appli", now(), sep="_")), ".csv"))
write.csv(cvs_f, file_name, row.names = FALSE)

cvs_m <- data.frame(cbind(Male_test$raw[1 : m_test], alord_m$cv, lord_m$cv, rhoalord_m$cv, rholord_m$cv))
names(cvs_m)[1] = "pvalues"
names(cvs_m)[2] = "alord"
names(cvs_m)[3] = "lord"
names(cvs_m)[4] = "rho_alord"
names(cvs_m)[5] = "rho_lord"
file_name = gsub(" " , "", paste("data/mfdr/impc/", gsub(" ", "_", paste("impc_male_mfdr_appli", now(), sep="_")), ".csv"))
write.csv(cvs_m, file_name, row.names = FALSE)

# --------------------------- Create plots and save them -------------------------------------------------------------
m_plot = param.list$m_plot

# -log(-log) function used to transform the y-scale in the plots 
llog <- function(x) { 
  -log(-log(x))
}
# female 

# nb_rej_alord <- sum(cvs_f$pvalues <= cvs_f$alord)
# nb_rej_lord <- sum(cvs_f$pvalues <= cvs_f$lord)
# nb_rej_rhoalord <- sum(cvs_f$pvalues <= cvs_f$rho_alord)
# nb_rej_rholord <- sum(cvs_f$pvalues <= cvs_f$rho_lord)
# nb_rej <- c(nb_rej_lord, nb_rej_rholord, nb_rej_alord, nb_rej_rhoalord)
# nb_rej

data <- data.frame(seq(1, m_plot), llog(cvs_f$pvalues[1 : m_plot]), llog(cvs_f$lord[1 : m_plot]), llog(cvs_f$rho_lord[1 : m_plot]))
names(data)[1] = "index"
names(data)[2] = "pvalues"
names(data)[3] = "lord_cv"
names(data)[4] = "rholord_cv"

data_rug_1 <- data.frame(seq(1, m_plot)[cvs_f$pvalues[1 : m_plot] <= cvs_f$rho_lord[1 : m_plot]])
names(data_rug_1)[1] <- "rholorddisco"

data_rug_2 <- data.frame(seq(1, m_plot)[cvs_f$pvalues[1 : m_plot] <= cvs_f$lord[1 : m_plot]])
names(data_rug_2)[1] <- "lorddisco"

plot_1_fem_loglog <- ggplot(data) +
  geom_point(aes_string("index", "pvalues"), size = 2, color = "black", alpha = 0.3) +
  geom_line(aes_string("index", "lord_cv"), lwd = 1, color = "darksalmon") +
  geom_line(aes_string("index", "rholord_cv"), lwd = 1, color = "darkolivegreen") +
  geom_rug(data = data_rug_1, mapping = aes(x = rholorddisco), color = "darkolivegreen", size = 1) +
  geom_rug(data = data_rug_2, mapping = aes(x = lorddisco), color = "darksalmon", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "t", y = "p-values / critical values") +
  ylim(c(min(data$pvalues), 0)) +
  theme(axis.text = element_text(size = 12))

ggsave("figures/application/mfdr/impc_app_female_lord_rholord.png", width = 16, height = 10, units = "cm")

data <- data.frame(seq(1, m_plot), llog(cvs_f$pvalues[1 : m_plot]), llog(cvs_f$alord[1 : m_plot]), llog(cvs_f$rho_alord[1 : m_plot]))
names(data)[1] = "index"
names(data)[2] = "pvalues"
names(data)[3] = "alord_cv"
names(data)[4] = "rhoalord_cv"

ata_rug_1 <- data.frame(seq(1, m_plot)[cvs_f$pvalues[1 : m_plot] <= cvs_f$rho_alord[1 : m_plot]])
names(data_rug_1)[1] <- "rhoalorddisco"

data_rug_2 <- data.frame(seq(1, m_plot)[cvs_f$pvalues[1 : m_plot] <= cvs_f$alord[1 : m_plot]])
names(data_rug_2)[1] <- "alorddisco"

plot_2_fem_loglog <- ggplot(data) +
  geom_point(aes_string("index", "pvalues"), size = 2, color = "black", alpha = 0.3) +
  geom_line(aes_string("index", "alord_cv"), lwd = 1, color = "darksalmon") +
  geom_line(aes_string("index", "rhoalord_cv"), lwd = 1, color = "darkolivegreen") +
  geom_rug(data = data_rug_1, mapping = aes(x = rhoalorddisco), color = "darkolivegreen", size = 1) +
  geom_rug(data = data_rug_2, mapping = aes(x = alorddisco), color = "darksalmon", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "t", y = "p-values / critical values") +
  ylim(c(min(data$pvalues), 0)) +
  theme(axis.text = element_text(size = 12))

ggsave("figures/application/mfdr/fem_alord_rhoalord.png", width = 16, height = 10, units = "cm")

# male 

# nb_rej_alord <- sum(cvs_m$pvalues <= cvs_m$alord)
# nb_rej_lord <- sum(cvs_m$pvalues <= cvs_m$lord)
# nb_rej_rhoalord <- sum(cvs_m$pvalues <= cvs_m$rho_alord)
# nb_rej_rholord <- sum(cvs_m$pvalues <= cvs_m$rho_lord)
# nb_rej <- c(nb_rej_lord, nb_rej_rholord, nb_rej_alord, nb_rej_rhoalord)
# nb_rej

data <- data.frame(seq(1, m_plot), llog(cvs_m$pvalues[1 : m_plot]), llog(cvs_m$lord[1 : m_plot]), llog(cvs_m$rho_lord[1 : m_plot]))
names(data)[1] = "index"
names(data)[2] = "pvalues"
names(data)[3] = "lord_cv"
names(data)[4] = "rholord_cv"

data_rug_1 <- data.frame(seq(1, m_plot)[cvs_m$pvalues[1 : m_plot] <= cvs_m$rho_lord[1 : m_plot]])
names(data_rug_1)[1] <- "rholorddisco"

data_rug_2 <- data.frame(seq(1, m_plot)[cvs_m$pvalues[1 : m_plot] <= cvs_m$lord[1 : m_plot]])
names(data_rug_2)[1] <- "lorddisco"

plot_1_male_loglog <- ggplot(data) +
  geom_point(aes_string("index", "pvalues"), size = 2, color = "black", alpha = 0.3) +
  geom_line(aes_string("index", "lord_cv"), lwd = 1, color = "darksalmon") +
  geom_line(aes_string("index", "rholord_cv"), lwd = 1, color = "darkolivegreen") +
  geom_rug(data = data_rug_1, mapping = aes(x = rholorddisco), color = "darkolivegreen", size = 1) +
  geom_rug(data = data_rug_2, mapping = aes(x = lorddisco), color = "darksalmon", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "t", y = "p-values / critical values") +
  ylim(c(min(data$pvalues), 0)) +
  theme(axis.text = element_text(size = 12))

ggsave("figures/application/mfdr/male_lord_rholord.png", width = 16, height = 10, units = "cm")

data <- data.frame(seq(1, m_plot), llog(cvs_m$pvalues[1 : m_plot]), llog(cvs_m$alord[1 : m_plot]), llog(cvs_m$rho_alord[1 : m_plot]))
names(data)[1] = "index"
names(data)[2] = "pvalues"
names(data)[3] = "alord_cv"
names(data)[4] = "rhoalord_cv"

data_rug_1 <- data.frame(seq(1, m_plot)[cvs_m$pvalues[1 : m_plot] <= cvs_m$rho_alord[1 : m_plot]])
names(data_rug_1)[1] <- "rhoalorddisco"

data_rug_2 <- data.frame(seq(1, m_plot)[cvs_m$pvalues[1 : m_plot] <= cvs_m$alord[1 : m_plot]])
names(data_rug_2)[1] <- "alorddisco"

plot_2_male_cvs_mfdr_loglog <- ggplot(data) +
  geom_point(aes_string("index", "pvalues"), size = 2, color = "black", alpha = 0.3) +
  geom_line(aes_string("index", "alord_cv"), lwd = 1, color = "darksalmon") +
  geom_line(aes_string("index", "rhoalord_cv"), lwd = 1, color = "darkolivegreen") +
  geom_rug(data = data_rug_1, mapping = aes(x = rhoalorddisco), color = "darkolivegreen", size = 1) +
  geom_rug(data = data_rug_2, mapping = aes(x = alorddisco), color = "darksalmon", outside = TRUE, size = 1) +
  coord_cartesian(clip = "off") +
  labs(x = "t", y = "p-values / critical values") +
  ylim(c(min(data$pvalues), 0)) +
  theme(axis.text = element_text(size = 12))

ggsave("figures/application/mfdr/male_alord_rhoalord.png", width = 16, height = 10, units = "cm")