#' kernel_bandwidth_mfdr
#'
#' Function used to compare the performance of rhoLORD and rhoALORD for
#' different rectangular kernel bandwidths.
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
#'                         Here the default is "rectangular" since we are studying the bandwidth for rectangular kernels.
#' @param q_2 An integer > 1 giving the bandwidth of the rectangular kernel.
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
#' @example kernel_bandwidth_mfdr(0.2, 0.1, 1000, 25, 100, 0.3, 0.4, "end", 0.5, "q-serie", 1.6, 10)
#'
#'
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
kernel_bandwidth_mfdr <- function(alpha, w0, nb_run, N, m, non_nulls_proportion,
                                          p3, cluster_option, lambda,
                                          gamma_type, q_1,
                                          kernel_bandwidth,
                                          gamma_prime_type = "rectangular",
                                          p1 = 0.01, p2 = 0.1) {

  gamma <- gamma_sequence(gamma_type, m, q_1)
  rectangular_gamma_prime <- gamma_sequence(gamma_prime_type, m, kernel_bandwidth)

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

                   # compute cv for rhoLORD and rhoALORD
                   proc_rho_lord <- rho_LORD(alpha, w0, raw.pvalues, CDF, gamma, rectangular_gamma_prime)
                   proc_rho_alord <- rho_ALORD(alpha, w0, raw.pvalues, CDF, gamma, lambda, rectangular_gamma_prime)

                   # estimate power and fwer for  rho_lord
                   output_2 <- number_of_discoveries(proc_rho_lord$rej, alternative_index, "mFDR")
                   rho_lord <- output_2$ratio_true_discoveries
                   fd_rho_lord <- output_2$error_quantity
                   nb_disco_rho_lord <- length(proc_rho_lord$rej)

                   # estimate power and fwer for  rho_alord
                   output_2 <- number_of_discoveries(proc_rho_alord$rej, alternative_index, "mFDR")
                   rho_alord <- output_2$ratio_true_discoveries
                   fd_rho_alord <- output_2$error_quantity
                   nb_disco_rho_alord <- length(proc_rho_alord$rej)

                   c(rho_lord, rho_alord,
                     fd_rho_lord, fd_rho_alord,
                     nb_disco_rho_lord, nb_disco_rho_alord)
                 }

  stopCluster(cl)
  colnames(res) <- NULL

  rho_lord <- data.frame(procedure = "rhoLORD",    ratio_true_discoveries = res[1, ])
  rho_alord <- data.frame(procedure = "rhoALORD",    ratio_true_discoveries = res[2, ])
  fd_rho_lord <- res[3, ]
  fd_rho_alord <- res[4, ]
  nb_disco_rho_lord <- pmax(res[5, ], 1)
  nb_disco_rho_alord <- pmax(res[6, ], 1)

  data_power <- data.frame(rbind(rho_alord, rho_lord))
  data_mfdr <- data.frame(procedure = c("rhoALORD", "rhoLORD"),
                          mFDR_estimate = c(mean(fd_rho_alord) / mean(nb_disco_rho_alord),
                                            mean(fd_rho_lord) / mean(nb_disco_rho_lord)
                          ))

  output <- list(data_power = data_power, data_mfdr = data_mfdr)
  return(output)
}

#' kernel_bandwidth_fwer
#'
#' Function used to compare the performance of rhoOB and rhoAOB for
#' different rectangular kernel bandwidths.
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
#' @param gamma_prime_type A string, either "log-q-serie", "q-serie", or "rectangular" to specify the type of the gamma prime sequence.
#'                         Here the default is "rectangular" since we are studying the bandwidth for rectangular kernels.
#' @param q_2 An integer > 1 giving the bandwidth of the rectangular kernel.
#' @param p1 A numeric corresponding to the Bernouilli parameter for generating a first group of nulls.
#'           Defaults to 0.01.
#' @param p2 A numeric corresponding to the Bernouilli parameter for generating a second group of nulls.
#'           Default to 0.1.
#' @param tau A numeric in [0, 1], used as a  threshold for hypotheses to be selected for testing.
#'            Used for Addis-spending only, kept Tian, J. and Ramdas, A. (2021) default values 0.5.
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
#' @example kernel_bandwidth_fwer(0.2, 1000, 25, 100, 0.3, 0.4, "end", 0.5, "q-serie", 1.6, 10)
#'
#' @import foreach
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster detectCores
#' @export
kernel_bandwidth_fwer <- function(alpha, nb_run, N, m, non_nulls_proportion,
                                              p3, cluster_option, lambda,
                                              gamma_type, q_1,
                                              kernel_bandwidth,
                                              gamma_prime_type = "rectangular",
                                              p1 = 0.01, p2 = 0.1) {

  gamma <- gamma_sequence(gamma_type, m, q_1)
  rectangular_gamma_prime <- gamma_sequence(gamma_prime_type, m, kernel_bandwidth)

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

                   # compute cv
                   proc_rhoob <- rho_OB(alpha, raw.pvalues, CDF, gamma, rectangular_gamma_prime)
                   proc_rhoaob <- rho_AOB(alpha, raw.pvalues, CDF, gamma, lambda, rectangular_gamma_prime)

                   # estimate power and fwer for  rhoob
                   output_2 <- number_of_discoveries(proc_rhoob$rej, alternative_index, "FWER")
                   rhoob <- output_2$ratio_true_discoveries
                   fd_rhoob <- output_2$error_quantity

                   # estimate power and fwer for  rhoaob
                   output_2 <- number_of_discoveries(proc_rhoaob$rej, alternative_index, "FWER")
                   rhoaob <- output_2$ratio_true_discoveries
                   fd_rhoaob <- output_2$error_quantity

                   c(rhoob, rhoaob, fd_rhoob, fd_rhoaob)
                 }

  stopCluster(cl)
  colnames(res) <- NULL

  rhoob <- data.frame(procedure = "rhoOB",    ratio_true_discoveries = res[1, ])
  rhoaob <- data.frame(procedure = "rhoAOB",    ratio_true_discoveries = res[2, ])
  fd_rhoob <- res[3, ]
  fd_rhoaob <- res[4, ]

  data_power <- data.frame(rbind(rhoaob, rhoob))
  data_fwer <- data.frame(procedure = c("rhoAOB", "rhoOB"),
                          FWER_estimate = c(mean(fd_rhoaob), mean(fd_rhoob)))

  output <- list(data_power = data_power, data_fwer = data_fwer)
  return(output)
}
