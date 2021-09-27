#include <Rcpp.h>
using namespace Rcpp;

//' stepfun
//'
//' Fast step function evaluation.
//' Step function is represented by a single numeric vector under the conditions
//' a) f(x) = x and b) it is complete, meaning that the last value is 1,
//' this is much faster than passing and evaluating R step function objects.
//'
// [[Rcpp::export]]
double stepfun(double &x, const NumericVector &sfun){
  // index variables and vector length
  int pos = 0, size = sfun.length();

  // computing results
  if(x < sfun[0]) return 0.0;
  if(x >= sfun[size - 1]) return sfun[size - 1];
  while(pos < size - 1 && x > sfun[pos]) pos++;
  if(sfun[pos] == x) return x;
  else return sfun[pos - 1];
}


//' rho_OB_AOB_Rcpp
//'
//' Function that can compute rho-OB and rho-AOB.
//' To compute rho_OB, lambda must be 0 (i.e. we don't threshold for adaptivity to the number of true nulls)
//' When gamma_prime is not provided and 'greedy' set to true, the greedy version of the procedure is computed.
//'
//' @param alpha  A numeric for the desired level of type I error control.
//' @param raw.pvalues A vector containing raw pvalues.
//' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
//'        Each support is represented by a vector in increasing order.
//' @param gamma A vector: the gamma spending sequence.
//' @param lambda A numeric in [0, 1] : the threshold for adaptivity.
//' @param gamma_prime A vector: the gamma prime smoothing sequence (it can be the same as gamma).
//' @param greedy A boolean to indicate whether we want to perform the greedy version (TRUE) or not (FALSE).
//'        Default to TRUE.
//'
//' @return A list containing a vector of the sequence of critical values and
//'         a vector of the indices of rejected hypothesis.
//'
//'
// [[Rcpp::export]]
List rho_OB_AOB_Rcpp(const double alpha, const NumericVector raw_pvalues, const List pCDFlist,
                     const NumericVector gamma, const double lambda,
                     const NumericVector gamma_prime = NumericVector(), const bool greedy = true)
{

  //-----------------------------------------------------------------------------
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1.");
  if(lambda < 0 || lambda > 1) stop("lambda must be between 0 and 1.");
  //-----------------------------------------------------------------------------

  int m = raw_pvalues.length();

  IntegerVector over_threshold(m);
  IntegerVector T(m);
  for(int i = 0; i < m; i++){
    over_threshold[i] = (int)(raw_pvalues[i] > lambda);
    if(i > 0) T[i] = T[i - 1] + over_threshold[i - 1];
  }

  // define critical values
  NumericVector cv = gamma[T];
  cv = cv * alpha * (1 - lambda) ;
  NumericVector u = clone(cv);

  for(int i = 1; i < m; i++){
    int end = i - 1;
    if (!over_threshold[i - 1]) {
        end = i - 2;
        cv[i] += cv[i - 1] - u[i - 1];
    }
    if(!greedy) {
        for(int j = 0; j <= end; j++){
        cv[i] += over_threshold[j] * gamma_prime[i - j - 1] * (cv[j] - stepfun(cv[j], as<NumericVector>(pCDFlist[j]))); // super uniformity reward
        }
    }else{
        cv[i] += over_threshold[i - 1]  * (cv[i - 1] - stepfun(cv[i - 1], as<NumericVector>(pCDFlist[i - 1]))); // only last super uniformity reward
    }
  }

  // rejection
  IntegerVector rej;
  for(int i = 0; i < m; i++) if(raw_pvalues[i] <= cv[i]) rej.push_back(i + 1);

  return List::create(Named("cv") = cv, Named("rej") = rej);
}


//' rho_LORD_ALORD_Rcpp
//'
//' Function that can compute rho-LORD and rho-ALORD.
//' To compute rho_LORD, lambda must be 0 (i.e. we don't threshold for adaptivity to the number of true nulls)
//' When gamma_prime is not provided and 'greedy' set to true, the greedy version of the procedure is computed.
//'
//' @param alpha  A numeric for the desired level of type I error control.
//' @param w0 A numeric in [0, alpha] representing the amount of earning back one wants.
//' @param raw.pvalues A vector containing raw pvalues.
//' @param pCDFlist A list containing the support of the discrete CDF of the p-values.
//'        Each support is represented by a vector in increasing order.
//' @param gamma A vector: the gamma spending sequence.
//' @param lambda A numeric in [0, 1] : the threshold for adaptivity.
//' @param gamma_prime A vector: the gamma prime smoothing sequence (it can be the same as gamma). Default to NULL.
//' @param greedy A boolean to indicate whether we want to perform the greedy version (TRUE) or not (FALSE).
//'        Default to TRUE.
//'
//' @return A list containing a vector of the sequence of critical values and
//'         a vector of the indices of rejected hypothesis.

//'
//'
// [[Rcpp::export]]
List rho_LORD_ALORD_Rcpp(const double alpha, const double w0, const NumericVector raw_pvalues,
                    const List pCDFlist, const NumericVector gamma, const double lambda,
                    const NumericVector gamma_prime = NumericVector(), const bool greedy = true)
{
  //-----------------------------------------------------------------------------
  if(alpha < 0 || alpha > 1) stop("alpha must be between 0 and 1.");
  if(lambda < 0 || lambda > 1) stop("lambda must be between 0 and 1.");
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
  NumericVector u = clone(cv);

  IntegerVector tau;
  if (raw_pvalues[0] <= cv[0]) tau.push_back(1);

  int T_1 = 0, T_j = 0;
  for(int i = 1; i < m; i++){
    double rej_reward_1 = 0, rej_reward_2 = 0;

    cv[i] += (cv[i - 1] - u[i - 1]) * (1 - over_threshold[i - 1]);
    if(greedy) {
      cv[i] += over_threshold[i - 1] * (cv[i - 1] - stepfun(cv[i - 1], as<NumericVector>(pCDFlist[i - 1])));
    }

    if(tau.length()){
      if(i == tau[0]) T_1 = 0;
      else if (i > tau[0]) T_1 = sum(over_threshold[Range(tau[0], i - 1)]);

      rej_reward_1 = (1 - lambda) * (alpha - w0) * gamma[T_1];
    }
    cv[i] += rej_reward_1;

    for(int j = 0; j <= i - 1; j++){
      if(!greedy) {
        cv[i] += over_threshold[j] * gamma_prime[i - j - 1] * (cv[j] - stepfun(cv[j], as<NumericVector>(pCDFlist[j]))); // super uniformity reward
      }


      if(j > 0 && j < tau.length()){
        if(i == tau[j]) T_j = 0;
        else if (i > tau[j]) T_j = sum(over_threshold[Range(tau[j], i - 1)]);

        rej_reward_2 += gamma[T_j];
      }
    }

    rej_reward_2 *= alpha * (1 - lambda);
    cv[i] += rej_reward_2;
    u[i] += rej_reward_1 + rej_reward_2;

    if(raw_pvalues[i] <= cv[i]) tau.push_back(i + 1);
  }

  // output
  return List::create(Named("cv") = cv, Named("rej") = tau);
}