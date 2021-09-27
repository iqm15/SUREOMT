#' data_simulation
#'
#' Function that generates data for performing Fisher's two-sided exact test for marginal counts.
#' Detailed description on the generation of data can be found in the reference given below.
#'
#' @param N An integer corresponding to the number of subjects studied (or the number of rows in the matrice).
#' @param m An integer corresponding to the number of hypotheses to test (or the number of columns in the matrice).
#' @param non_nulls_proportion A numeric in [0, 1] corresponding to the quantity of signal the user wants in the data.
#' @param p3 A numeric in [0, 1] corresponding to the strength of the signal the user wants.
#' @param cluster_option Either "end", "begin", "begin_middle", "begin_end", "midlle_end", or "no_cluster_shuffle".
#'        This option indicates how to position the signal in the stream of hypothesis.
#' @param p1 A numeric corresponding to the Bernouilli parameter for generating a first group of nulls.
#'           Default to 0.01.
#' @param p2 A numeric corresponding to the Bernouilli parameter for generating a second group of nulls.
#'           Default to 0.1.
#' @param seed An integer for setting the seed to reproduce exactly a simulation.
#'
#' @return A list containing the simulated data as a data frame and the indices of the non-nulls in the data frame.
#'         The data frame has 4 columns of length m. The first and third columns represent
#'         group_1 and group_2, each of their component represents a count in {0...N}.
#'         The second and fourth columns represent the total number N which is the same for
#'         group_1 and group_2.
#'
#' @examples data_simulation(25, 100, 0.3, 0.4, "end").
#'
#' @references Sebastian Döhler, Guillermo Durand, Etienne Roquain
#'             "New FDR bounds for discrete and heterogeneous tests,"
#'             \emph{Electronic Journal of Statistics, Electron. J. Statist. 12(1), 1867-1900, (2018).}
#'             \url{https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-12/issue-1/New-FDR-bounds-for-discrete-and-heterogeneous-tests/10.1214/18-EJS1441.full}
#'
#' @export
data_simulation <- function(N, m, non_nulls_proportion, p3,
                            cluster_option, p1 = 0.01, p2 = 0.1, seed=NULL) {

  #-----------------------------------------------------------------------------
  if (!is.null(seed)){
    set.seed(seed)
  }
  #-----------------------------------------------------------------------------
  testthat::expect_gt(N, 0)
  testthat::expect_gt(m, 0)
  testthat::expect_lte(p1, 1)
  testthat::expect_lte(p2, 1)
  testthat::expect_lte(p3, 1)
  testthat::expect_lte(non_nulls_proportion, 1)
  match.arg(cluster_option, c("end", "begin", "begin_end", "begin_middle",
                              "middle_end", "no_cluster_shuffle"))


  proportions <- c((1 - non_nulls_proportion) / 2,
                   (1 - non_nulls_proportion) / 2,
                   non_nulls_proportion)
  quantities <- round(m * proportions)
  try (if (sum(round(quantities)) != m) stop("the proportion you provided give unrounded quantities,
                                             you should change the proportion"))
  #-----------------------------------------------------------------------------

  probabilities_1 <- c(p1, p2, p2)
  probabilities_2 <- c(p1, p2, p3)


  if (cluster_option == "end"){
    new_quantities <- quantities
    probs_1 <- probabilities_1
    probs_2 <- probabilities_2
    alternative_index <- seq(m - quantities[3] + 1, m)
  }

  if (cluster_option == "begin") {
    new_quantities <- c(quantities[3], quantities[2], quantities[1])
    probs_1 <- c(probabilities_1[3], probabilities_1[2], probabilities_1[1])
    probs_2 <- c(probabilities_2[3], probabilities_2[2], probabilities_2[1])
    alternative_index <- seq(1, new_quantities[1])
  }

  if (cluster_option == "begin_end") {
    half = quantities[3] / 2
    new_quantities <- c(half, quantities[1], quantities[2], half)
    probs_1 <- c(probabilities_1[3], probabilities_1[1], probabilities_1[2], probabilities_1[3])
    probs_2 <- c(probabilities_2[3], probabilities_2[1], probabilities_2[2], probabilities_2[3])

    alternative_index <- c(seq(1, half), seq(m - half + 1, m))
  }

  if (cluster_option == "begin_middle") {
    half = quantities[3] / 2
    new_quantities <- c(half, quantities[1], half, quantities[2])
    probs_1 <- c(probabilities_1[3], probabilities_1[1], probabilities_1[3], probabilities_1[2])
    probs_2 <- c(probabilities_2[3], probabilities_2[1], probabilities_2[3], probabilities_2[2])

    alternative_index <- c(seq(1, half), seq(sum(new_quantities[1:2]) + 1, sum(new_quantities[1:3])))
  }

  if (cluster_option == "middle_end") {
    half = quantities[3] / 2
    new_quantities <- c(quantities[1], half, quantities[2], half)
    probs_1 <- c(probabilities_1[1], probabilities_1[3], probabilities_1[2], probabilities_1[3])
    probs_2 <- c(probabilities_2[1], probabilities_2[3], probabilities_2[2], probabilities_2[3])

    alternative_index <- c(seq(new_quantities[1] + 1, sum(new_quantities[1:2])),
                           seq(sum(new_quantities[1:3]) + 1, m))
  }

  if (cluster_option == "no_cluster_shuffle") {
    new_quantities <- quantities
    probs_1 <- probabilities_1
    probs_2 <- probabilities_2
  }

  data_matrix1 <- c()
  data_matrix2 <- c()

  for (i in 1:length(new_quantities)){
    size = new_quantities[i] * N
    data_matrix1 <- c(data_matrix1, rbinom(size, 1, probs_1[i]))
    data_matrix2 <- c(data_matrix2, rbinom(size, 1, probs_2[i]))
  }

  group_1 <- matrix(data = data_matrix1, nrow = m, ncol = N, byrow = TRUE)
  group_2 <- matrix(data = data_matrix2, nrow = m, ncol = N, byrow = TRUE)

  x1 <- c(rowSums(group_1))
  x2 <- c(rowSums(group_2))
  y1 <- rep(N, m)
  y2 <- rep(N, m)

  if (cluster_option == "no_cluster_shuffle") {
    shuffled <- shuffle_vec(x1)
    x1 <- shuffled$shuffle_vec
    x2 <- shuffle_vec(x2, shuffled$permutation_index)$shuffle_vec
    alternative_index = shuffled$permutation_index[(m - quantities[3] + 1): m]
  }

  # verify that we have the right nb of alternative index
  testthat::expect_length(alternative_index, quantities[3])

  data <- data.frame(x1, y1, x2, y2)
  output <- list(data = data, alternative_index = alternative_index)
  return(output)
}


#' pvalues-simulation
#'
#' Function that performs Fisher's exact two-sided test.
#' This function uses the 'fisher.pvalues.support' from
#' the DiscreteFDR package of Junge. F et al (2019).
#'
#' @param data frame A data frame of 4 columns and any number of rows (the number of rows
#'                  is the number of hypotheses that the user is testing).
#' @param alternative A string informing if we are doing two.sided test or not.
#'                    The argument can either be "two.sided", "less", or "greater".
#'                    Here the default is to "two.sided".
#' @param input The format of the input data frame, see the references below.
#'              Here the default is to "marginal".
#'
#' @return A list containing the raw p-values and their CDF support.
#'
#' @examples test_data <- data.frame(c(1, 0, 2, 5, 2), rep(25, 5), c(0, 0, 2, 14, 11), rep(25, 5))
#'           or test_data <- data_simulation(25, 100, 0.3, 0.4, "end")$data
#'           pvalues_simulation(test_data).
#'
#' @references Junge. F, Durand. G, Döhler. S, and Roquain. E (2019)
#'             DiscreteFDR: An R package for controlling the false
#'             discovery rate for discrete test statistics.
#'             \url{https://arxiv.org/pdf/1904.02054.pdf}
#'             \url{https://cran.r-project.org/web/packages/DiscreteFDR/index.html}.
#'
#'
#' @export
pvalues_simulation <- function(dataframe, alternative = "two.sided", input = "marginal") {

  #-----------------------------------------------------------------------------
  if (!is.data.frame(dataframe)) {
    stop("You must provide a dataframe to compute the p-values.")
  }
  #-----------------------------------------------------------------------------

  p <- DiscreteFDR::fisher.pvalues.support(dataframe, alternative, input)

  return(p)
}
