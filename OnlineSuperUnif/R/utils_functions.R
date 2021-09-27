#'
#'
#'@export
get_normalizingconstant <- function(type, q, N = 1000) {
  #-----------------------------------------------------------------------------
  match.arg(type, c("log-q-serie", "q-serie", "JM"))
  #-----------------------------------------------------------------------------
  i <- 1:N
  if (type == "q-serie") {
    gamma <- 1 / ((i)^q)
    integral_upper_bound = sum(gamma) - ((1 / (1 - q)) * (N)^(1-q))
  }

  else if (type == "log-q-serie"){
    gamma <- 1 / ((i + 1) * (log(i + 1)^q))
    integral_upper_bound = sum(gamma) - (1 / ((1 - q) * (log(N + 1)^(1 - q))))
  }

  else{
    gamma <- log(pmax(i, 2)) / (i * exp(sqrt(log(i))))
    integral_upper_bound = sum(gamma) + (2 * exp(-sqrt(log(N))) * (log(N)^(3 / 2) + 3 * log(N) + 6 * sqrt(log(N)) + 6))
  }

  return(1 / integral_upper_bound)

}


#' gamma_sequence.
#'
#' Function that computes a nonnegative decreasing sequence.
#' The user can choose to make the sequence sum to exactly one
#' (and thus using the number of hypotheses to test),
#' or to make the sequence sum to less than one by approximating the infinity.
#' Three choices for the type of sequence are proposed, of which
#' log-q serie and q-serie as proposed by Tian and Ramdas (2021).
#'
#' @param type Either "log-q-serie", "q-serie" or a "rectangular" kernel.
#' @param nb_pvalues An integer giving the nb of p-values (/ hypothesis) to test.
#' @param q The exponent for computing the sequence or the kernel bandwidth.
#'          Note that when using a rectangular kernel, q must be an integer.
#'
#' @return A vector: the gamma sequence.
#'
#' @example gamma_sequence("log-q-serie", 100, 2).
#'
#' @references Tian, J. and Ramdas, A. (2021). Online control of the familywise
#'             error rate. \emph{Statistical Methods for Medical Research},
#'                          \url{https://journals.sagepub.com/eprint/AYRRKZX7XMTVHKCFYBJY/full}
#'
#' @export
gamma_sequence <- function(type, nb_pvalues, q) {

  #-----------------------------------------------------------------------------
  match.arg(type, c("log-q-serie", "q-serie", "JM", "rectangular"))
  if (type == "rectangular"){
    if (q %% 1 != 0) {
      stop("For using a rectangular kernel, you should provide an integer for the bandwidth q")
    }
  }
  #-----------------------------------------------------------------------------
  if (type != "rectangular"){
    normalization_constant = get_normalizingconstant(type, q)
  }

  if (type == "log-q-serie") {
    i <- 1:nb_pvalues
    gamma <- 1 / ((i + 1) * (log(i + 1)^q))

    # normalize the sequence
    gamma = gamma * normalization_constant
  }

  else if (type == "q-serie") {
    i <- 1:nb_pvalues
    gamma <- 1 / ((i)^q)

    # normalize the sequence
    gamma = gamma * normalization_constant

  }

  else if (type == "JM") {
    i <- 1:nb_pvalues
    gamma <- log(pmax(i, 2)) / (i * exp(sqrt(log(i))))

    # normalize the sequence
    gamma = gamma * normalization_constant
  }

  else {
    if (q - round(q) != 0) {
      stop("you should provide a round number for the bandwidth, q, when wantingto use a rectangular kernel")
    }
    gamma <- c(rep(1 / q, q), rep(0, nb_pvalues - q))
  }

  testthat::expect_lte(sum(gamma), 1) # test that the sum is less than or equal to 1
  return(gamma)

}


#' shuffle_vec
#'
#' Function that shuffles a vector (permutation).
#' This function allows to study the signal position scheme where
#' the signal is not clustered but positioned randomly across the whole stream of hypothesis
#' (signal position = "no_cluster_shuffle" in data_simulation function).
#'
#'
#' @param vec A vector that needs to be shuffled.
#' @param permutation_index A vector indicating how to shuffle the vector
#'                          if one wants to perform a certain permutation.
#'
#' @return A list containing the shuffled vector and the index of the entries.
#'
#' @example shuffle_vec(c(11, 12, 13, 14, 15), c(4, 3, 5, 1, 2))
#'          should return the permuted vector c(14, 13, 15, 11, 12),
#'           and c(4, 3, 5, 1, 2), the permutation index.
#'
#' @export
shuffle_vec <- function(vec, permutation_index = NULL) {


  if (missing(permutation_index)){
    l = length(vec)
    permutation_index <- gtools::permute(1:l)
  }

  permutation_mat <- as.matrix(Matrix::sparseMatrix(seq_along(permutation_index),
                                            permutation_index, x=1))
  shuffled_vec <- as.vector(vec %*% permutation_mat)

  output <- list(shuffle_vec = shuffled_vec, permutation_index = permutation_index)
  return(output)
}


#' number_of_discoveries
#'
#' Function that allows to get the necessary quantities to estimate the
#' error (power, FWER or mFDR).
#'
#' @param rej_index A vector containing the indices of the rejected hypothesis.
#' @param alternative_index A vector containing the indices (in the stream of hypothesis) of the signal.
#' @param error_metric A string, either "FWER" or "mFDR" to indicate the error metric the user is studying.
#'
#' @return A list containing
#'         ratio_true_discoveries : Ratio between the nb of true discoveries
#'                                  and the number of non-nulls (= signals),
#'         Nb of true discoveries,
#'         error_quantity : depending on the error metric;
#'                          either a boolean stating the presence of a false discovery (FWER),
#'                          or the number of false discoveries (mFDR).
#'
#' @example number_of_discoveries(c(4, 5, 13, 14, 17), seq(13, 20), "FWER") should
#'          return (3 / 20, 3, 1) (where 1 stands for TRUE) and
#'          number_of_discoveries(c(4, 5, 13, 14, 17), seq(13, 20), "mFDR") should
#'          return (3 / 20, 3, 2)
#'
#' @export
number_of_discoveries <- function(rej_index, alternative_index, error_metric) {
  #-----------------------------------------------------------------------------
  match.arg(error_metric, c("mFDR", "FWER"))
  #-----------------------------------------------------------------------------

  nb_true_discoveries <- sum(rej_index %in% alternative_index)

  ratio_true_discoveries <- nb_true_discoveries / length(alternative_index)

  if (error_metric == "FWER") {
    false_discoveries_bool <- (length(rej_index) > nb_true_discoveries)
    error_quantity <- false_discoveries_bool
  }

  else if (error_metric == "mFDR") {
    nb_false_discoveries <- length(rej_index) - nb_true_discoveries
    error_quantity <- nb_false_discoveries
  }

  output <- list(ratio_true_discoveries = ratio_true_discoveries,
                 nb_true_discoveries = nb_true_discoveries,
                 error_quantity = error_quantity)

  return(output)

}


#' get_CDF
#'
#' Function that allows getting the CDF of p-values ready to plot.
#' This function is used only for shiny apps.
#'
#' @param N An integer corresponding to the number of subjects studied (or the number of rows in the matrice).
#' @param m An integer corresponding to the number of hypotheses to test (or the number of columns in the matrice).
#' @param non_nulls_proportion A numeric in [0, 1] corresponding to the quantity of signal the user wants in the data.
#' @param p3 A numeric in [0, 1] corresponding to the strength of the signal the user wants.
#' @param cluster_option Either "end", "begin", "begin_middle", "begin_end", "midlle_end", or "no_cluster_shuffle".
#'                       This option indicates how to position the signal in the stream of hypothesis.
#' @param p1 A numeric corresponding to the Bernouilli parameter for generating a first group of nulls.
#' @param p2 A numeric corresponding to the Bernouilli parameter for generating a second group of nulls.
#'
#' @return A list with the p-values' CDFs ready to use.
#'
#' @example get_CDF(25, 100, 0.3, 0.4, "end").
#'
get_CDF <- function(N, m, non_nulls_proportion, p3,
                    cluster_option, p1 = 0.01, p2 = 0.1) {

  proportions = c((1 - non_nulls_proportion) / 2,
                  (1 - non_nulls_proportion) / 2,
                  non_nulls_proportion)

  data <- data_simulation(N, m, non_nulls_proportion, p3, cluster_option)$data

  CDF_list <- pvalues_simulation(data)$support
  stepf <- lapply(CDF_list, function(x) stepfun(x, c(0, x)))
  return(stepf)

}


#' male_female_pvalue_min
#'
#'
male_female_pvalue_min <- function(Male_test, Fem_test) {
  pvalues <- numeric(nrow(Male_test))
  support <- list(nrow(Male_test))

  for (i in 1:nrow(Male_df)) {
    if (min(Male_test$raw[i], Fem_test$raw[i]) == Male_test$raw[i]) {
      pvalues[i] = Male_test$raw[i]
      support[i] = Male_test$support[i]
    }
    else {
      pvalues[i] = Fem_test$raw[i]
      support[i] = Fem_test$support[i]
    }
  }
  output <- list(raw = pvalues, support = support)
}


#' get_counts
#'
#' @export
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


