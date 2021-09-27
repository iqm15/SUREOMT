#' addis_onlineFDR
#'
#' Function that can compute the ADDIS procedure.
#' This function used the function ADDIS of
#' the OnlineFDR package, see reference below.
#' For the parameters set with default values we used
#' the values given by Tian. J, and Ramdas. A (2019).
#'
#' @param alpha  A numeric for the desired level of type I error control.
#' @param w0 A numeric representing the initial `wealth' of the procedure (the default values is \eqn{\tau \lambda
#'          \alpha/2}
#' @param raw.pvalues A vector containing the raw p-values.
#' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
#'                 Each support is represented by a vector in increasing order.
#' @param gamma A vector: the gamma spending sequence.
#' @param lambda A numeric in [0, 1]: the threshold for adaptivity.
#'               Default to 0.25.
#' @param tau: A numeric in [0, 1]: the threshold for hypotheses to be selected for testing.
#'              Default to 0.5.
#'
#' @return A list containing a vector of the sequence of critical values and
#'         a vector of the indices of rejected hypothesis.
#'
#' @example test_data <- data_simulation(25, 100, 0.3, 0.4, "end")
#'          test <- pvalues_simulation(test_data$data)
#'          raw.pvalues <- test$raw
#'          gamma <- gamma_sequence("q-serie", 100, 1.6,)
#'          addis_onlineFDR(0.2, raw.pvalues, gamma)
#'
#'
#' @references Tian, J. and Ramdas, A. (2019). ADDIS: an adaptive discarding
#'             algorithm for online FDR control with conservative nulls.
#'             \emph{Advances in Neural Information Processing Systems}, 9388-9396.
#' @references Robertson DS, Liou L, Ramdas A, Javanmard A, Tian J, Zrnic T, Karp NA (2019).
#'             onlineFDR: Online error control. R package 2.1.0.
#'             \url{https://bioconductor.org/packages/devel/bioc/html/onlineFDR.html}
#'
#'
#' @export
addis_onlineFDR <- function(alpha, w0, raw.pvalues, gamma, lambda = 0.25, tau = 0.5) {
  out <- onlineFDR::ADDIS(raw.pvalues, alpha, gamma, w0, lambda, tau)
  R <- out$R
  cv <- out$alphai
  rej <- which(out$R == 1)
  output <- list(cv = cv, rej = rej)

  return(output)
}


#' addis_spending_onlineFDR
#'
#' Function that can compute ADDIS-spending procedure.
#' This function used the function ADDIS_spending of
#' the OnlineFDR package, see reference below.
#' For the parameters set with default values we used
#' the values given by Tian. J, and Ramdas. A (2021).
#'
#' @param alpha  A numeric for the desired level of type I error control.
#' @param raw.pvalues A vector containing the raw p-values.
#' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
#'                 Each support is represented by a vector in increasing order.
#' @param gamma A vector: the gamma spending sequence.
#' @param lambda A numeric in [0, 1]: the threshold for adaptivity.
#'               Default to 0.25.
#' @param tau: A numeric in [0, 1]: the threshold for hypotheses to be selected for testing.
#'              Default to 0.5.
#'
#' @return A list containing a vector of the sequence of critical values and
#'         a vector of the indices of rejected hypothesis.
#'
#' @example test_data <- data_simulation(25, 100, 0.3, 0.4, "end")
#'          test <- pvalues_simulation(test_data$data)
#'          raw.pvalues <- test$raw
#'          gamma <- gamma_sequence("q-serie", 100, 1.6,)
#'          addis_spending_onlineFDR(0.2, raw.pvalues, gamma)
#'
#'
#' @references Tian, J. and Ramdas, A. (2021). Online control of the familywise
#'             error rate. \emph{Statistical Methods for Medical Research},
#'                          \url{https://journals.sagepub.com/eprint/AYRRKZX7XMTVHKCFYBJY/full}.
#' @references Robertson DS, Liou L, Ramdas A, Javanmard A, Tian J, Zrnic T, Karp NA (2019).
#'             onlineFDR: Online error control. R package 2.1.0.
#'             \url{https://bioconductor.org/packages/devel/bioc/html/onlineFDR.html}
#'
#'
#' @export
addis_spending_onlineFDR <- function(alpha, raw.pvalues, gamma, lambda = 0.25, tau = 0.5) {
  out <- onlineFDR::ADDIS_spending(raw.pvalues, alpha, gamma, lambda, tau)
  R <- out$R
  cv <- out$alphai
  rej <- which(out$R == 1)
  output <- list(cv = cv, rej = rej)

  return(output)
}


#' OB
#'
#' Function that can compute OB procedure.
#'
#' @param alpha  A numeric for the desired level of type I error control.
#' @param raw.pvalues A vector containing the raw p-values.
#' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
#'                 Each support is represented by a vector in increasing order.
#' @param gamma A vector: the gamma spending sequence.
#'
#' @return A list containing a vector of the sequence of critical values and
#'         a vector of the indices of rejected hypothesis.
#'
#' @example test_data <- data_simulation(25, 100, 0.3, 0.4, "end")
#'          test <- pvalues_simulation(test_data$data)
#'          raw.pvalues <- test$raw
#'          gamma <- gamma_sequence("q-serie", 100, 1.6,)
#'          OB(0.2, raw.pvalues, gamma)
#'
#' @export
OB <- function(alpha, raw.pvalues, gamma) {
  return(OB_AOB_Rcpp(alpha, raw.pvalues, gamma, lambda = 0))
}


#' AOB
#'
#' Function that can compute AOB (Adaptive-Spending of Tian & Ramdas (2021)) procedure.
#'
#' @param alpha  A numeric for the desired level of type I error control.
#' @param raw.pvalues A vector containing the raw p-values.
#' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
#'                 Each support is represented by a vector in increasing order.
#' @param gamma A vector: the gamma spending sequence.
#' @param lambda A numeric in [0, 1] : the adaptivity threshold.
#'
#' @return A list containing a vector of the sequence of critical values and
#'         a vector of the indices of rejected hypothesis.
#'
#' @example test_data <- data_simulation(25, 100, 0.3, 0.4, "end")
#'          test <- pvalues_simulation(test_data$data)
#'          raw.pvalues <- test$raw
#'          gamma <- gamma_sequence("q-serie", 100, 1.6,)
#'          lambda = 0.5
#'          AOB(0.2, raw.pvalues, gamma, lambda)
#'
#'
#'
#' @export
AOB <- function(alpha, raw.pvalues, gamma, lambda) {
  return(OB_AOB_Rcpp(alpha, raw.pvalues, gamma, lambda))
}


#' rho_OB
#'
#' Function that can compute rho-OB procedure.
#'
#' @param alpha  A numeric for the desired level of type I error control.
#' @param raw.pvalues A vector containing the raw p-values.
#' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
#'                 Each support is represented by a vector in increasing order.
#' @param gamma A vector: the gamma spending sequence.
#' @param gamma_prime A vector: the gamma prime smoothing sequence (it can be the same as gamma).
#'                    When gamma_prime is not provided, the greedy version of the procedure is performed.
#'
#' @return A list containing a vector of the sequence of critical values and
#'         a vector of the indices of rejected hypothesis.
#'
#' @example test_data <- data_simulation(25, 100, 0.3, 0.4, "end")
#'          test <- pvalues_simulation(test_data$data)
#'          raw.pvalues <- test$raw
#'          CDF <- test$support
#'          gamma <- gamma_sequence("q-serie", 100, 1.6,)
#'          gamma_prime = gamma
#'          rho_OB(0.2, raw.pvalues, CDF, gamma, gamma_prime)
#'
#' @export
rho_OB <- function(alpha, raw.pvalues, pCDFlist, gamma, gamma_prime) {
  if (missing(gamma_prime)){
    warning("The gamma prime sequence is not provided,
             the greedy version of the procedure will be computed")
    return(rho_OB_AOB_Rcpp(alpha, raw.pvalues, pCDFlist, gamma, lambda = 0))
  } else {
    return(rho_OB_AOB_Rcpp(alpha, raw.pvalues, pCDFlist, gamma, lambda = 0, gamma_prime, greedy = FALSE))
  }
}


#' rho_AOB
#'
#' Function that can compute rho-AOB procedure.
#'
#' @param alpha  A numeric for the desired level of type I error control.
#' @param raw.pvalues A vector containing the raw p-values.
#' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
#'                 Each support is represented by a vector in increasing order.
#' @param gamma A vector: the gamma spending sequence.
#' @param gamma_prime A vector: the gamma prime smoothing sequence (it can be the same as gamma).
#'                    When gamma_prime is not provided, the greedy version of the procedure is performed.
#' @param lambda A numeric in [0, 1] : the threshold for adaptivity.
#'
#' @return A list containing a vector of the sequence of critical values and
#'         a vector of the indices of rejected hypothesis.
#'
#' @example test_data <- data_simulation(25, 100, 0.3, 0.4, "end")
#'          test <- pvalues_simulation(test_data$data)
#'          raw.pvalues <- test$raw
#'          CDF <- test$support
#'          gamma <- gamma_sequence("q-serie", 100, 1.6,)
#'          gamma_prime = gamma
#'          lambda = 0.5
#'          rho_AOB(0.2, raw.pvalues, CDF, gamma, lamnda, gamma_prime)
#'
#' @references Tian, J. and Ramdas, A. (2021). Online control of the familywise
#'             error rate. \emph{Statistical Methods for Medical Research},
#'                          \url{https://journals.sagepub.com/eprint/AYRRKZX7XMTVHKCFYBJY/full}.
#'
#' @export
rho_AOB <- function(alpha, raw.pvalues, pCDFlist, gamma, lambda,  gamma_prime) {
  if (missing(gamma_prime)) {
    warning("The gamma prime sequence is not provided,
             the greedy version of the procedure will be computed")
    return(rho_OB_AOB_Rcpp(alpha, raw.pvalues, pCDFlist, gamma, lambda))
  } else {
    return(rho_OB_AOB_Rcpp(alpha, raw.pvalues, pCDFlist, gamma, lambda, gamma_prime, greedy = FALSE))
  }
}


#' rho_LORD
#'
#' Function that can compute rho-LORD procedure.
#'
#' @param alpha  A numeric in [0, 1] for the desired level of type I error control.
#' @param w0 A numeric representing the initial wealth.
#' @param raw.pvalues A vector containing the raw p-values.
#' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
#'                 Each support is represented by a vector in increasing order.
#' @param gamma A vector: the gamma spending sequence.
#' @param gamma_prime A vector: the gamma prime smoothing sequence (it can be the same as gamma).
#'                    When gamma_prime is not provided, the greedy version of the procedure is performed.
#'
#' @return A list containing a vector of the sequence of critical values and
#'         a vector of the indices of rejected hypothesis.
#'
#' @example test_data <- data_simulation(25, 100, 0.3, 0.4, "end")
#'          test <- pvalues_simulation(test_data$data)
#'          raw.pvalues <- test$raw
#'          CDF <- test$support
#'          gamma <- gamma_sequence("q-serie", 100, 1.6,)
#'          gamma_prime = gamma
#'          rho_LORD(0.2, raw.pvalues, CDF, gamma, lamnda, gamma_prime)
#'
#'
#' @export
rho_LORD <- function(alpha, w0, raw.pvalues, pCDFlist, gamma, gamma_prime) {
  if (missing(gamma_prime)) {
    warning("The gamma prime sequence is not provided,
             the greedy version of the procedure will be computed")
    return(rho_LORD_ALORD_Rcpp(alpha, w0, raw.pvalues, pCDFlist, gamma, lambda = 0))
  } else {
    return(rho_LORD_ALORD_Rcpp(alpha, w0, raw.pvalues, pCDFlist, gamma, lambda = 0, gamma_prime, greedy = FALSE))
  }
}


#' rho_ALORD
#'
#' Function that can compute rho-ALORD procedure.
#'
#' @param alpha  A numeric in [0, 1] for the desired level of type I error control.
#' @param w0 A numeric representing the initial wealth.
#' @param raw.pvalues A vector containing the raw p-values.
#' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
#'                 Each support is represented by a vector in increasing order.
#' @param gamma A vector: the gamma spending sequence.
#' @param gamma_prime A vector: the gamma prime smoothing sequence (it can be the same as gamma).
#'                    When gamma_prime is not provided, the greedy version of the procedure is performed.
#'
#' @return A list containing a vector of the sequence of critical values and
#'         a vector of the indices of rejected hypothesis.
#'
#' @example test_data <- data_simulation(25, 100, 0.3, 0.4, "end")
#'          test <- pvalues_simulation(test_data$data)
#'          raw.pvalues <- test$raw
#'          CDF <- test$support
#'          gamma <- gamma_sequence("q-serie", 100, 1.6,)
#'          gamma_prime = gamma
#'          rho_ALORD(0.2, raw.pvalues, CDF, gamma, lamnda, gamma_prime)
#'
#'
#' @export
rho_ALORD <- function(alpha, w0, raw.pvalues, pCDFlist, gamma, lambda, gamma_prime){
  if (missing(gamma_prime)) {
    warning("The gamma prime sequence is not provided,
             the greedy version of the procedure will be computed")
    return(rho_LORD_ALORD_Rcpp(alpha, w0, raw.pvalues, pCDFlist, gamma, lambda))
  } else {
    return(rho_LORD_ALORD_Rcpp(alpha, w0, raw.pvalues, pCDFlist, gamma, lambda, gamma_prime, greedy = FALSE))
  }
}
