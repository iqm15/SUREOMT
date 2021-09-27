#' compare_fwer_proc_on_impc
#'
#' @export
compare_fwer_proc_on_impc <- function(alpha, nb_run, m,
                                      lambda,
                                      gamma_type, q_1,
                                      gamma_prime_type, q_2) {

  # compute the gamma and gamma_prime sequence
  print(q_1)
  print(q_2)

  gamma <- gamma_sequence(gamma_type, m, q_1)
  gamma_prime <- gamma_sequence(gamma_prime_type, m, q_2)

  print("done gamma")
  # load the data sets
  Male_df <- read.csv("/users/home/meah/phd_projects/online_superunif_mt/code/my_impc_data/impc_female_df.csv")[1 : m, ]
  Fem_df <- read.csv("/users/home/meah/phd_projects/online_superunif_mt/code/my_impc_data/impc_female_df.csv")[1 : m, ]

  print("done")
  # create CPU cluster for parallel computing
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)
  print("cluster creation done")
  # for reproducible random numbers on parallel processors
  # registerDoRNG(seed = 123)

  res <- foreach(i = 1:nb_run, .combine = cbind,
                 .export = c("male_female_pvalue_min"),
                 .packages = c("testthat", "DiscreteFDR"))%dopar%{
                  print("in foreach loop")
                  Male_test <- fisher.pvalues.support(Male_df[sample(1:m, m, replace = FALSE), 1:ncol(Male_df)],
                                                      alternative = "greater", input = "noassoc")
                  Fem_test <- fisher.pvalues.support(Fem_df[sample(1:m, m, replace = FALSE), 1:ncol(Fem_df)],
                                                     alternative = "greater", input = "noassoc")
                  print("getting pvalues done")

                  # for each gene select the lowest p-value between the one of female and the one of male
                  selection <- male_female_pvalue_min(Male_test, Fem_test)
                  pvalues <- selection$raw
                  CDF <- selection$support

                  # perform procedures
                  OB_ <- OB(alpha, pvalues, gamma)
                  AOB_ <- AOB(alpha, pvalues, gamma, lambda)
                  rho_OB_ <- rho_OB(alpha, pvalues, CDF, gamma, gamma_prime)
                  rho_AOB_ <- rho_AOB(alpha, pvalues, CDF, gamma, lambda, gamma_prime)

                  # collect the number of discoveries for each procedures
                  nb_disco_AOB = length(AOB_$rej)
                  nb_disco_OB = length(OB_$rej)
                  nb_disco_rho_AOB = length(rho_AOB_$rej)
                  nb_disco_rho_OB = length(rho_OB_$rej)

                  c(nb_disco_AOB, nb_disco_OB, nb_disco_rho_AOB, nb_disco_rho_OB)
                 }

  stopCluster(cl)
  colnames(res) <- NULL

  aob  <- data.frame(procedure = "AOB", discoveries = res[1, ])
  ob  <- data.frame(procedure = "OB", discoveries = res[2, ])
  rho_aob  <- data.frame(procedure = "RHO_AOB", discoveries = res[3, ])
  rho_ob  <- data.frame(procedure = "RHO_OB", discoveries = res[4, ])

  data_discoveries <- data.frame(rbind(aob, ob, rho_aob, rho_ob))

  return()
}
