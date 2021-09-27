#include <Rcpp.h>
using namespace Rcpp;

//' OB_AOB_Rcpp
//'
//' Function that can compute OB and AOB (Adaptive-Spending of Tian & Ramdas (2021)) procedure.
//' Setting lambda to 0 gives OB.
//'
//' @param alpha A numeric in [0, 1] for the desired level of type I error control.
//' @param raw.pvalues A vector containing raw p-values.
//' @param gamma A vector : the gamma spending sequence.
//' @param lambda A numeric in [0, 1] : the adaptivity threshold that sets the
//'          threshold for 'candidate' hypothesis. Default to 0.5,
//'          setting lambda to 0 gives 0B the non adaptive versio of AOB.
//'
//' @references Tian, J. and Ramdas, A. (2021). Online control of the familywise
//'             error rate. \emph{Statistical Methods for Medical Research},
//'                          \url{https://journals.sagepub.com/eprint/AYRRKZX7XMTVHKCFYBJY/full}.
//'
//'
// [[Rcpp::export]]
List OB_AOB_Rcpp(const double alpha, const NumericVector raw_pvalues,
            const NumericVector gamma, const double lambda) {
  //-----------------------------------------------------------------------------
  if (alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1");
  if (lambda < 0 || lambda > 1) stop("lambda must be in [0, 1]");
  //-----------------------------------------------------------------------------

  int m = raw_pvalues.length();
  LogicalVector over_threshold(m);
  IntegerVector T(m);
  for(int i = 0; i < m; i++){
    over_threshold[i] = (raw_pvalues[i] > lambda);
    if(i > 0){
      T[i] = T[i - 1] + (int)over_threshold[i - 1];
    }
  }

  // critical values
  NumericVector cv(m);
  for(int i = 0; i < m; i++){
    cv[i] = alpha * (1 - lambda) * gamma[T[i]];
  }

  // rejections
  IntegerVector rej;
  for(int i = 0; i < m; i++) if(raw_pvalues[i] <= cv[i]) rej.push_back(i + 1);

  return List::create(Named("cv") = cv, Named("rej") = rej);

}


//' lord_OnlineSupeUnif
//'
//' Function that computes LORD of Ramdas et al. (2017).
//'
//' @param alpha  A numeric for the desired level of type I error control.
//' @param w0 A numeric in [0, alpha] representing the amount of earning back one wants.
//' @param raw.pvalues A vector containing raw pvalues.
//' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
//'        Each support is represented by a vector in increasing order.
//' @param gamma A vector: the gamma spending sequence.
//'
//' @return A list containing a vector of the sequence of critical values and
//'         a vector of the indices of rejected hypothesis.
//'
//'
//' @references Ramdas, A., Yang, F., Wainwright M.J. and Jordan, M.I. (2017). Online control
//'             of the false discovery rate with decaying memory. \emph{Advances in Neural
//'             Information Processing Systems 30}, 5650-5659.
//'
//'
//' @export
// [[Rcpp::export]]
List lord_OnlineSuperUnif (const double alpha, const double w0,
                           const NumericVector raw_pvalues, const NumericVector gamma)
{
  //-----------------------------------------------------------------------------
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1.");
  if (w0 > alpha) stop("The wealth must be less or equal to alpha");
  //-----------------------------------------------------------------------------

  int m = raw_pvalues.length();

  // define critical values
  NumericVector cv = clone(gamma);
  cv = cv * w0;

  IntegerVector tau;
  int k = 0;
  if (raw_pvalues[0] <= cv[0]){
    tau.push_back(1);
    k++;
  }

  int T_1 = 0, T_j = 0;
  for(int i = 1; i < m; i++){
    double rej_reward_1 = 0, rej_reward_2 = 0;

    if(k > 0)
      rej_reward_1 = (alpha - w0) * gamma[i - tau[0]];

    cv[i] += rej_reward_1;

    for(int j = 1; j < k; j++)
      rej_reward_2 += gamma[i - tau[j]];

    cv[i] += alpha * rej_reward_2;

    if(raw_pvalues[i] <= cv[i]){
      tau.push_back(i + 1);
      k++;
    }
  }

  // output
  return List::create(Named("cv") = cv, Named("rej") = tau);
}


//' saffron_OnlineSuperUnif
//'
//' Function that computes SAFFRON of Ramdas et al. (2018).
//'
//' @param alpha  A numeric for the desired level of type I error control.
//' @param w0 A numeric in [0, alpha] representing the amount of earning back one wants.
//' @param raw.pvalues A vector containing raw p-values.
//' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
//'        Each support is represented by a vector in increasing order.
//' @param gamma A vector: the gamma spending sequence.
//' @param gamma_prime A vector: the gamma prime smoothing sequence (it can be the same as gamma). Default to NULL.
//' @param lambda A numeric in [0, 1] : the threshold for adaptivity.
//'
//' @return A list containing a vector of the sequence of critical values and
//'         a vector of the indices of rejected hypothesis.
//'
//' @reference Ramdas, A., Zrnic, T., Wainwright M.J. and Jordan, M.I. (2018).
//'             SAFFRON: an adaptive algorithm for online control of the false discovery
//'             rate. \emph{Proceedings of the 35th International Conference in Machine
//'             Learning}, 80:4286-4294.
//'
//'
//' @export
// [[Rcpp::export]]
List saffron_OnlineSuperUnif(const double alpha, const double w0, const NumericVector raw_pvalues,
                              const NumericVector gamma, const double lambda, const bool capping = false)
{
  //-----------------------------------------------------------------------------
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1.");
  if(lambda < 0 || lambda > 1) stop("lambda must be between 0 and 1.");
  if (w0 > alpha * (1- lambda)) stop("The wealth must be less or equal to alpha * lambda");
  //-----------------------------------------------------------------------------

  int m = raw_pvalues.length();

  IntegerVector over_threshold(m);
  IntegerVector T_0(m);
  for(int i = 0; i < m; i++){
    over_threshold[i] = (int)(raw_pvalues[i] > lambda);
    if(i > 0) T_0[i] = T_0[i - 1] + over_threshold[i - 1];
  }

  // define critical values
  NumericVector cv = gamma[T_0];
  cv = cv * (1 - lambda) * w0;
  if (capping) cv[0] = std::min(cv[0], lambda);


  IntegerVector tau;
  if (raw_pvalues[0] <= cv[0]) tau.push_back(1);

  double rej_reward_1, rej_reward_2;
  int T_1 = 0, T_j = 0;
  for(int i = 1; i < m; i++){
    rej_reward_1 = 0;
    rej_reward_2 = 0;

    if(tau.length()){
      if(i == tau[0]) T_1 = 0;
      else if (i > tau[0]) T_1 = sum(over_threshold[Range(tau[0], i - 1)]);

      rej_reward_1 = (1 - lambda) * (alpha - w0) * gamma[T_1];
    }
    cv[i] += rej_reward_1;

    for(int j = 1; j <= i - 1; j++){

      if(j > 0 && j < tau.length()){
        if(i == tau[j]) T_j = 0;
        else if (i > tau[j]) T_j = sum(over_threshold[Range(tau[j], i - 1)]);

        rej_reward_2 += gamma[T_j];
      }
    }

    rej_reward_2 *= alpha * (1 - lambda);
    cv[i] += rej_reward_2;

    if(capping){
      cv[i] = std::min(lambda, cv[i]);
    }

    if(raw_pvalues[i] <= cv[i]) tau.push_back(i + 1);
  }

  // output
  return List::create(Named("cv") = cv, Named("rej") = tau);
}
