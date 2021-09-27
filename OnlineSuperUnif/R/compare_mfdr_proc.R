#' compare_mfdr_proc
#'
#' Function that collects data to estimate the power and mFDR of procedures
#' LORD, ALORD, SAFFRON, rho-LORD, rho-ALORD, and ADDIS to compare their performance.
#'
#' @param alpha A numeric between 0 and 1 giving the overall error budget.
#' @param w0 A numeric between 0 and alpha, the initial wealth of the procedures.
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
#' @param gamma_prime_type A string, either "log-q-serie", "q-serie", or "rectangular" to specify the type of the gamma prime sequence.
#' @param q_2 A numeric > 1 indicating the exponent for computing the gamma prime sequence.
#'              Note that for a rectangular kernel, q_2 must be an integer.
#' @param p1 A numeric corresponding to the Bernouilli parameter for generating a first group of nulls.
#'           Defaults to 0.01.
#' @param p2 A numeric corresponding to the Bernouilli parameter for generating a second group of nulls.
#'           Default to 0.1.
#' @param tau A numeric in [0, 1], used as a threshold for hypotheses to be selected for testing.
#'            Used for ADDIS only, kept Tian, J. and Ramdas, A. (2019) default values 0.5.
#'
#' @return Returns two data frames.
#'         - The first data frame contains the ratio between the number of true discoveries and the number of non-nulls
#'           for each run performed for each procedure.
#'           This ratio must be averaged over all runs, for each procedure to estimate the power.
#'         - The second data frame contains an estimate of the mFDR for each procedure i.e. the ratio between
#'           the mean (over nb_run) of the number of false discoveries
#'           and the mean (over nb_run) of the number of discoveries.
#'
#'
#' @example compare_mfdr_proc(0.2, 0.1, 1000, 25, 100, 0.3, 0.4, "end", 0.5, "q-serie", 1.6, "log-q-serie", 2)
#'
#'
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
compare_mfdr_proc <- function(alpha, w0, nb_run, N, m, non_nulls_proportion,
                              p3, cluster_option, lambda,
                              gamma_type, q_1,
                              gamma_prime_type, q_2,
                              p1 = 0.01, p2 = 0.1) {

  gamma <- gamma_sequence(gamma_type, m, q_1)
  gamma_prime <- gamma_sequence(gamma_prime_type, m, q_2)

  # create CPU cluster for parallel computing
  cl <- makeCluster(parallel::detectCores() - 1)
  registerDoParallel(cl)

  # for reproducible random numbers on parallel processors
  #registerDoRNG(seed = 1)

  res <- foreach(i = 1:nb_run, .combine = cbind,
                 .export = c("data_simulation", "pvalues_simulation", "number_of_discoveries", "shuffle_vec"),
                 .packages = c("testthat", "DiscreteFDR", "gtools", "Matrix"))%dopar%{

                   # generate data
                   output_1 <- data_simulation(N, m, non_nulls_proportion, p3, cluster_option)
                   data <- output_1$data
                   alternative_index <- output_1$alternative_index

                   # run test and store p-values
                   test <- pvalues_simulation(data)
                   raw.pvalues <- test$raw
                   CDF <- test$support

                   # compute cv for LORD, ALORD, SAFFRON, rhoLORD, rhoALORD and ADDIS
                   proc_rho_lord <- rho_LORD(alpha, w0, raw.pvalues, CDF, gamma, gamma_prime)
                   proc_rho_alord <- rho_ALORD(alpha, w0, raw.pvalues, CDF, gamma, lambda, gamma_prime)

                   proc_lord <- lord_OnlineSuperUnif(alpha, w0, raw.pvalues, gamma)
                   proc_saffron <- saffron_OnlineSuperUnif(alpha, w0, raw.pvalues, gamma, lambda, capping = TRUE)
                   proc_alord <- saffron_OnlineSuperUnif(alpha, w0, raw.pvalues, gamma, lambda)

                   # for addis the value for tau and lambda are set to the default given by Tian and Ramdas
                   proc_addis <- addis_onlineFDR(alpha, (alpha * 0.25 * 0.5) / 2 , raw.pvalues, gamma, lambda = 0.25, tau = 0.5)

                   # estimate power and fwer for rho_lord
                   output_2 <- number_of_discoveries(proc_rho_lord$rej, alternative_index, "mFDR")
                   rho_lord <- output_2$ratio_true_discoveries
                   fd_rho_lord <- output_2$error_quantity
                   nb_disco_rho_lord <- length(proc_rho_lord$rej)

                   # estimate power and fwer for rho_alord
                   output_2 <- number_of_discoveries(proc_rho_alord$rej, alternative_index, "mFDR")
                   rho_alord <- output_2$ratio_true_discoveries
                   fd_rho_alord <- output_2$error_quantity
                   nb_disco_rho_alord <- length(proc_rho_alord$rej)

                   # estimate power and fwer for lord
                   output_2 <- number_of_discoveries(proc_lord$rej, alternative_index, "mFDR")
                   lord <- output_2$ratio_true_discoveries
                   fd_lord <- output_2$error_quantity
                   nb_disco_lord <- length(proc_lord$rej)

                   # estimate power and fwer for saffron
                   output_2 <- number_of_discoveries(proc_saffron$rej, alternative_index, "mFDR")
                   saffron <- output_2$ratio_true_discoveries
                   fd_saffron <- output_2$error_quantity
                   nb_disco_saffron <- length(proc_saffron$rej)

                   # estimate power and fwer for alord
                   output_2 <- number_of_discoveries(proc_alord$rej, alternative_index, "mFDR")
                   alord <- output_2$ratio_true_discoveries
                   fd_alord <- output_2$error_quantity
                   nb_disco_alord <- length(proc_alord$rej)

                   # estimate power and fwer for addis
                   output_2 <- number_of_discoveries(proc_addis$rej, alternative_index, "mFDR")
                   addis <- output_2$ratio_true_discoveries
                   fd_addis <- output_2$error_quantity
                   nb_disco_addis <- length(proc_addis$rej)

                   c(rho_lord, rho_alord, lord, saffron, alord, addis,
                     fd_rho_lord, fd_rho_alord, fd_lord, fd_saffron, fd_alord, fd_addis,
                     nb_disco_rho_lord, nb_disco_rho_alord, nb_disco_lord, nb_disco_saffron, nb_disco_alord, nb_disco_addis)
                 }
  stopCluster(cl)
  colnames(res) <- NULL

  rho_lord    <- data.frame(procedure = "rhoLORD",    ratio_true_discoveries = res[1, ])
  rho_alord <- data.frame(procedure = "rhoALORD", ratio_true_discoveries = res[2, ])
  lord     <- data.frame(procedure = "LORD",     ratio_true_discoveries = res[3, ])
  saffron  <- data.frame(procedure = "SAFFRON",  ratio_true_discoveries = res[4, ])
  alord <- data.frame(procedure = "ALORD",  ratio_true_discoveries = res[5, ])
  addis <- data.frame(procedure = "ADDIS",  ratio_true_discoveries = res[6, ])
  fd_rho_lord <- res[7, ]
  fd_rho_alord <- res[8, ]
  fd_lord     <- res[9, ]
  fd_saffron  <- res[10, ]
  fd_alord <- res[11, ]
  fd_addis <- res[12, ]
  nb_disco_rho_lord    <- pmax(res[13, ], 1)
  nb_disco_rho_alord <- pmax(res[14, ], 1)
  nb_disco_lord     <- pmax(res[15, ], 1)
  nb_disco_saffron  <- pmax(res[16, ], 1)
  nb_disco_alord <- pmax(res[17, ], 1)
  nb_disco_addis <- pmax(res[18, ], 1)

  data_power <- data.frame(rbind(addis, alord, lord, rho_alord, rho_lord, saffron))
  data_mfdr <- data.frame(procedure = c("ADDIS", "ALORD", "LORD", "rhoALORD", "rhoLORD", "SAFFRON"),
                          mFDR_estimate = c(mean(fd_addis) / mean(nb_disco_addis),
                                            mean(fd_alord) / mean(nb_disco_alord),
                                            mean(fd_lord) / mean(nb_disco_lord),
                                            mean(fd_rho_alord) / mean(nb_disco_rho_alord),
                                            mean(fd_rho_lord) / mean(nb_disco_rho_lord),
                                            mean(fd_saffron) / mean(nb_disco_saffron)
                          ))

  output <- list(data_power = data_power, data_mfdr = data_mfdr)
  return(output)
}
