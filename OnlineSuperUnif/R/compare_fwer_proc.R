#' compare_fwer_proc
#'
#' Function that collects data to estimate the power and FWER of procedure
#' OB, AOB, rho-OB, rho-AOB, and Addis-spending to compare their performance.
#'
#' @param alpha A numeric between 0 and 1 giving the overall error budget.
#' @param nb_run An Integer giving the number of runs to make. If we put 100 runs, the
#'               estimated power and fwer will be done upon a mean of 100 quantities.
#' @param N An integer corresponding to the number of subjects studied (or the number of rows in the matrice).
#' @param m An integer corresponding to the number of hypotheses to test (or the number of columns in the matrice).
#' @param non_nulls_proportion A numeric in [0, 1] corresponding to the quantity of signal the user wants in the data.
#' @param p3 A numeric in [0, 1] corresponding to the strength of the signal the user wants.
#' @param cluster_option Either "end", "begin", "begin_middle", "begin_end", "midlle_end", or "no_cluster_shuffle".
#'                       This option indicates how to position the signal in the stream of hypothesis.
#' @param lambda A numeric in [0, 1], to threshold the p-value for the adaptive procedures.
#' @param gamma_type A string, either "log-q-serie", "q-serie" or "rectangular" to specify the type of the gamma sequence.
#' @param q_1 A numeric > 1 indicating the exponent for computing the gamma sequence.
#' @param gamma_prime_type A string, either "log-q-serie", "q-serie" or "rectangular" to specify the type of the gamma prime sequence.
#' @param q_2 A numeric > 1 indicating the exponent for computing the gamma prime sequence.
#'            Note that when using a rectangular kernel, q_2 must be an integer.
#' @param p1 A numeric corresponding to the Bernouilli parameter for generating a first group of nulls.
#'           Default to 0.01.
#' @param p2 A numeric corresponding to the Bernouilli parameter for generating a second group of nulls.
#'           Default to 0.1.
#' @param tau A numeric in [0, 1], used as a  threshold for hypotheses to be selected for testing.
#'            Used for Addis-spending only, kept Tian, J. and Ramdas, A. (2021) default values 0.5.
#'
#' @return Returns two data frames.
#'         - The first data frame contains the ratio between the number of true discoveries and the number of non-nulls
#'           for each run performed, for each procedure.
#'           This ratio must be averaged over all runs for each procedure to estimate the power.
#'         - The second data frame contains averaged (over nb_run) booleans, indicating whether
#'           a false discovery was made or not, for each procedure.
#'           It's the FWER estimate for each procedure.
#'
#'
#' @example compare_fwer_proc(0.2, 1000, 25, 100, 0.3, 0.4, "end", 0.5, "q-serie", 1.6, "log-q-serie", 2)
#'
#'
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
compare_fwer_proc <- function(alpha, nb_run, N, m, non_nulls_proportion,
                              p3, cluster_option, lambda,
                              gamma_type, q_1,
                              gamma_prime_type,
                              q_2,
                              p1 = 0.01, p2 = 0.1) {

  gamma <- gamma_sequence(gamma_type, m, q_1)
  gamma_prime <- gamma_sequence(gamma_prime_type, m, q_2)


  # create CPU cluster for parallel computing
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)

  # for reproducible random numbers on parallel processors
  #registerDoRNG(seed = 123)

  res <- foreach(i = 1:nb_run, .combine = cbind,
                 .export = c("data_simulation", "pvalues_simulation",
                             "number_of_discoveries", "shuffle_vec"), # is that necessary now that these are all in the package ?
                 .packages = c("testthat", "DiscreteFDR", "gtools", "Matrix"))%dopar%{

                   # generate data
                   output_1 <- data_simulation(N, m, non_nulls_proportion, p3, cluster_option)
                   data <- output_1$data
                   alternative_index <- output_1$alternative_index

                   # run test and store p-values
                   test <- pvalues_simulation(data)
                   raw.pvalues <- test$raw
                   CDF <- test$support

                   # perform procedures OB, AOB, rhoOB, rhoAOB, and ADDIS-spending
                   proc_rho_ob <- rho_OB(alpha, raw.pvalues, CDF, gamma, gamma_prime)
                   proc_rho_aob <- rho_AOB(alpha, raw.pvalues, CDF, gamma, lambda, gamma_prime)

                   proc_ob <- OB(alpha, raw.pvalues, gamma)
                   proc_aob <- AOB(alpha, raw.pvalues, gamma, lambda)

                   # For addis-spending the value for lambda and tau are hardcoded to the default values given by Tian and Ramdas
                   proc_addis_spend <- addis_spending_onlineFDR(alpha, raw.pvalues, gamma, lambda = 0.25, tau = 0.5)

                   # estimate power and fwer for rho_ob
                   output_2 <- number_of_discoveries(proc_rho_ob$rej, alternative_index, "FWER")
                   rho_ob <- output_2$ratio_true_discoveries
                   fd_rho_ob <- output_2$error_quantity

                   # estimate power and fwer for rho_aob
                   output_2 <- number_of_discoveries(proc_rho_aob$rej, alternative_index, "FWER")
                   rho_aob <- output_2$ratio_true_discoveries
                   fd_rho_aob <- output_2$error_quantity

                   # estimate power and fwer for aob
                   output_2 <- number_of_discoveries(proc_aob$rej, alternative_index, "FWER")
                   aob <- output_2$ratio_true_discoveries
                   fd_aob <- output_2$error_quantity

                   # estimate power and fwer for ob
                   output_2 <- number_of_discoveries(proc_ob$rej, alternative_index, "FWER")
                   ob <- output_2$ratio_true_discoveries
                   fd_ob <- output_2$error_quantity

                   # estimate power and fwer for addis-spending
                   output_2 <- number_of_discoveries(proc_addis_spend$rej, alternative_index, "FWER")
                   addis_spend <- output_2$ratio_true_discoveries
                   fd_addis_spend <- output_2$error_quantity

                   c(rho_ob, rho_aob, aob, ob, addis_spend,
                     fd_rho_ob, fd_rho_aob, fd_aob, fd_ob, fd_addis_spend)
                 }

  stopCluster(cl)
  colnames(res) <- NULL
  rho_ob  <- data.frame(procedure = "RHO_OB", ratio_true_discoveries = res[1, ])
  rho_aob  <- data.frame(procedure = "RHO_AOB", ratio_true_discoveries = res[2, ])
  aob  <- data.frame(procedure = "AOB", ratio_true_discoveries = res[3, ])
  ob  <- data.frame(procedure = "OB", ratio_true_discoveries = res[4, ])
  addis_spend <- data.frame(procedure = "ADDIS-spending", ratio_true_discoveries = res[5, ])
  fd_rho_ob   <- res[6, ]
  fd_rho_aob  <- res[7, ]
  fd_aob  <- res[8, ]
  fd_ob <- res[9, ]
  fd_addis_spend  <- res[10, ]

  data_power <- data.frame(rbind(addis_spend, aob, ob, rho_aob, rho_ob))
  data_fwer  <- data.frame(procedure = c("ADDIS-spending", "AOB", "OB", "rhoAOB", "rhoOB"),
                           fwer_estimate = c(mean(fd_addis_spend), mean(fd_aob),
                                             mean(fd_ob), mean(fd_rho_aob), mean(fd_rho_ob))
  )

  output <- list(data_power = data_power, data_fwer = data_fwer)
  return(output)
}
