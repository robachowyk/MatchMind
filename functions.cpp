#include <random>
#include <chrono>

#define ARMA_64BIT_WORD 1 
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]] 
using namespace Rcpp;

static std::random_device rd; 

// initialize Mersennes' twister using rd to generate the seed
static std::mt19937 gen{rd()}; 
std::uniform_real_distribution<double> dist(0, 1);

// [[Rcpp::export]]
NumericVector callFunction(NumericVector x, Function f) {
  NumericVector res = f(x);
  return res;
}

// [[Rcpp::export]] 
List F11(const List & F, const int & nvals)
{
  List out(nvals);
  for (int k = 0; k < nvals; k++)
  {
    List tmp = F[k];
    IntegerVector tmpA = tmp[0];
    IntegerVector tmpB = tmp[1];
    IntegerMatrix tmpC(tmpA.length() * tmpB.length(),2);
    int counter = 0;
    for (int i = 0; i < tmpA.length(); i++)
    {
      for (int j = 0; j < tmpB.length(); j++)
      {
        tmpC(counter,0) = tmpA[i];
        tmpC(counter,1) = tmpB[j];
        counter += 1;
      }
    }
    out[k] = tmpC;
  }
  return out;
}

// [[Rcpp::export]] 
List F1(const IntegerVector & HA, const IntegerVector & HB, const int & nvals)
{
  List out(nvals);
  for (int i = 0; i < nvals; i++)
  { IntegerVector tmpA;
    IntegerVector tmpB;
    List both(2);
    both[0] = tmpA;
    both[1] = tmpB;
    out[i] = both;
  }
  for (int i = 0; i < HA.length(); i++)
  {
    List tmp = out[HA[i]-1];
    IntegerVector tmpA = tmp[0];
    tmpA.push_back(i+1);
    tmp[0] = tmpA;
    out[HA[i]-1] = tmp;
  } 
  for (int j = 0; j < HB.length(); j++)
  {
    List tmp = out[HB[j]-1];
    IntegerVector tmpB = tmp[1];
    tmpB.push_back(j+1);
    tmp[1] = tmpB;
    out[HB[j]-1] = tmp;
  }
  return out;
}

// [[Rcpp::export]] 
List sampleD(const IntegerMatrix & S,
             const NumericVector & LLA, 
             const NumericVector & LLB, 
             const arma::sp_mat & LLL, 
             const NumericVector & gamma, 
             double loglik, 
             arma::sp_mat D, 
             int nlinkrec,
             LogicalVector sumRowD,
             LogicalVector sumColD)
{
  for (int q = 0; q < S.nrow(); q++)
  {
    int i = S(q,0)-1;
    int j = S(q,1)-1;
    // If non match and possible match -> check if match
    if((sumRowD(i)==false) && (sumColD(j)==false) && LLL(i,j) < 10000000)
    {					
      double loglikNew = loglik  
      // Comparison vectors
      - LLB(j) - LLA(i) 
      + LLL(i,j)
      // Bipartite matching
      - log(1-gamma(i)) + log(gamma(i))
      - log(LLB.length() - nlinkrec);
      double sumlogdensity = log(1 + exp(loglik-loglikNew)) + loglikNew;
      double pswitch = exp(loglikNew - sumlogdensity);		
      // Random number smaller than prob -> generate binomial value
      bool link = dist(gen) < pswitch;
      if(link)
      {
        loglik = loglikNew;
        D(i,j) = true;
        sumRowD(i) = true;
        sumColD(j) = true;
        nlinkrec = nlinkrec + 1; 
      }
    }else if(D(i,j)==true)
    {
      // If match -> check if nonmatch
      double loglikNew = loglik
      // Comparison vectors
      + LLB(j) + LLA(i) 
      - LLL(i,j)
      // Bipartite matching
      + log(1-gamma(i)) - log(gamma(i))
      + log(LLB.length() - nlinkrec+1);
      double sumlogdensity = log(1 + exp(loglik-loglikNew)) + loglikNew;	
      double pswitch = exp(loglikNew - sumlogdensity);		
      bool nolink = dist(gen) < pswitch;
      if(nolink)
      {
        loglik = loglikNew;
        D(i,j) = false;
        sumRowD(i) = false;
        sumColD(j) = false;
        nlinkrec = nlinkrec - 1; 
      }			
    }
  }
  arma::uvec indices = find(D);
  int n = indices.n_elem;
  IntegerMatrix links(n, 2);
  for (int i = 0; i < n; i++) {
    links(i, 0) = indices(i) % D.n_rows;
    links(i, 1) = indices(i) / D.n_rows;
  }
  // Return to R
  List ret;
  ret["D"] = D;
  ret["links"] = links;
  ret["sumRowD"] = sumRowD;	
  ret["sumColD"] = sumColD;
  ret["loglik"] = loglik;
  ret["nlinkrec"] = nlinkrec;
  return ret;
}

// [[Rcpp::export]]
IntegerVector sampleNL(IntegerVector G, NumericVector eta, NumericVector phi)
{
  IntegerVector H(G.length());
  // Number of possible values
  int nval = eta.length(); 
  // Create the possible values to sample from
  IntegerVector choice_set = seq_len(nval);
  // Possible registration errors
  double pMissing = phi[1];
  double pTypo = (1-pMissing) * (1-phi[0]) / (eta.length()-1);
  double pAgree = (1-pMissing) * phi[0];
  // Iterate over all elements
  for(int i = 0; i < G.length(); i++) 
  {
    // Create a vector indicating P(Registered=X|True)
    // First value is for the missings
    // What happens if missing:
    if(G(i) == 0)
    { 
      // All equally likely
      H(i) = Rcpp::sample(choice_set, 1, false, eta)[0];
    }else
    {
      NumericVector help1(nval, pTypo);   
      help1(G(i)-1) = pAgree; 
      // Joint probability to have the registered and true value
      NumericVector prob = eta * help1;  
      H(i) = Rcpp::sample(choice_set, 1, false, prob)[0];
    }
  }
  return H;
}

// [[Rcpp::export]]
IntegerVector sampleL(IntegerVector GA, IntegerVector GB, NumericVector survivalpSameH,
                      IntegerMatrix choice_set, IntegerVector choice_equal,
                      int nval, NumericVector phikA, NumericVector phikB, NumericVector eta)
{
  IntegerVector H(GA.length());
  int size_choice_set = choice_set.nrow();
  IntegerVector choice_index = seq_len(size_choice_set);
  // Possible actions for file A
  double pMissingA = phikA[1];
  double pTypoA = (1-pMissingA) * (1-phikA[0]) / (nval-1);
  double pAgreeA = (1-pMissingA) * phikA[0];
  // Possible actions for file B
  double pMissingB = phikB[1];
  double pTypoB = (1-pMissingB) * (1-phikB[0]) / (nval-1);
  double pAgreeB = (1-pMissingB) * phikB[0];
  // Iterate over all matches
  for(int i = 0; i < GA.length(); i++) 
  {
    // Prob that both TRUE values are the same
    double pSameH = survivalpSameH[i];
    // Dummy Vectors
    // Define P(Hb|Ha)
    NumericVector probH(size_choice_set, pSameH);
    //Define P(Ga|Ha) and P(Gb|Hb)
    NumericVector helpA(size_choice_set, pTypoA);   
    NumericVector helpB(size_choice_set, pTypoB);
    for(int j = 0; j < size_choice_set; j++) 
    {
      // Prob to observe HB=b|HA=a where a!=b
      if(choice_equal(j)==0)
      {
        probH(j) = (1-pSameH)/(nval-1); 
      }
      if(GA(i) == choice_set(j,0))
        helpA(j) = pAgreeA;
      if(GB(i) == choice_set(j,1))
        helpB(j) = pAgreeB;
    }
    // There are four options possible 
    // Both missing
    if(GA(i)==0 && GB(i)==0)
    {
      NumericVector prob = eta * probH;
      H(i) = Rcpp::sample(choice_index, 1, false, prob)[0];
    }else if(GA(i)>0 && GB(i)==0)
    {
      // Joint probability to have the registered and true value
      NumericVector prob = eta * probH * helpA;  
      H(i) = Rcpp::sample(choice_index, 1, false, prob)[0];
    }else if(GA(i)==0 && GB(i)>0)
    {
      // Joint probability to have the registered and true value
      NumericVector prob = eta * probH  * helpB;  
      H(i) = Rcpp::sample(choice_index, 1, false, prob)[0];
    }else if(GA(i)>0 && GB(i)>0)
    {  
      // None missing
      // Create vectors indicating P(Registered=X|True)
      // Joint probability to have the registered and true value
      NumericVector prob = eta * probH * helpA * helpB; 
      H(i) = Rcpp::sample(choice_index, 1, false, prob)[0];
    }   
  }
  return H;
}

// [[Rcpp::export]]
IntegerMatrix cartesianProduct(IntegerVector vec1, IntegerVector vec2) {
  int n = vec1.size();
  int m = vec2.size();
  IntegerMatrix result(n * m, 2);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m; j++) {
      result(i * m + j, 0) = vec1[i];
      result(i * m + j, 1) = vec2[j];
    }
  }
  return result;
}

// [[Rcpp::export]]
IntegerMatrix ExpandGrid(IntegerVector vector1, IntegerVector vector2) {
  return cartesianProduct(vector1, vector2);
}

// [[Rcpp::export]]
IntegerVector generateSequence(int n) {
  IntegerVector result(n);
  for (int i = 0; i < n; i++) {
    result[i] = i + 1;
  }
  return result;
}

// [[Rcpp::export]]
List sampleH(IntegerVector nA, IntegerVector nB, IntegerMatrix links, List omegaData, LogicalVector pivs_stable, List pivsA, List pivsB, IntegerVector nvalues, arma::sp_mat D, LogicalVector nonlinkedA, LogicalVector nonlinkedB, List eta, List omega, List phi)
{
  IntegerMatrix truepivsA(nA[0],nA[1]);
  IntegerMatrix truepivsB(nB[0],nB[1]);
  int nphi = 2;
  for (int k = 0; k < nvalues.length(); k++) {
    IntegerVector truepivsA_k = truepivsA(_,k);
    IntegerVector truepivsB_k = truepivsB(_,k);
    NumericVector eta_k = eta[k];
    NumericVector phi_k = phi[k];
    NumericVector phi_k_A(nphi);
    phi_k_A[0] = phi_k[0];
    phi_k_A[1] = phi_k[1];
    NumericVector phi_k_B(nphi);
    phi_k_B[0] = phi_k[0];
    phi_k_B[1] = phi_k[2];
    IntegerVector pivsA_k = pivsA[k];
    IntegerVector pivsB_k = pivsB[k];
    IntegerVector pivsA_k_L = pivsA_k[links(_,0)];
    IntegerVector pivsB_k_L = pivsB_k[links(_,1)];
    IntegerVector pivsA_k_NL = pivsA_k[nonlinkedA];
    IntegerVector pivsB_k_NL = pivsB_k[nonlinkedB];
    IntegerVector truepivsA_k_NL = sampleNL(pivsA_k_NL, eta_k, phi_k_A);
    IntegerVector truepivsB_k_NL = sampleNL(pivsB_k_NL, eta_k, phi_k_B);
    truepivsA_k[nonlinkedA] = truepivsA_k_NL;
    truepivsB_k[nonlinkedB] = truepivsB_k_NL;
    IntegerMatrix choice_set;
    IntegerVector choice_equal;
    NumericVector eta_choice;
    if (links.nrow()>0)
    {
      IntegerVector values = generateSequence(nvalues[k]);
      choice_set = ExpandGrid(values, values);
      choice_equal = choice_set(_,1) == choice_set(_,0);
      eta_choice = eta_k[choice_set(_,1) - 1];
      NumericVector survivalpSameH(links.nrow(), 1.0);
      if(!pivs_stable[k])
      {
        NumericVector omegaData_k = omegaData[k];
        survivalpSameH = exp(-omegaData_k);
      }
      IntegerVector out = sampleL(pivsA_k_L, pivsB_k_L, survivalpSameH, choice_set, choice_equal, nvalues[k], phi_k_A, phi_k_B, eta_choice);
      IntegerVector choice_set_A = choice_set(_,0);
      IntegerVector choice_set_B = choice_set(_,1);
      truepivsA_k[links(_,0)] = choice_set_A[out - 1];
      truepivsB_k[links(_,1)] = choice_set_B[out - 1];
    }
    truepivsA(_,k) = truepivsA_k;
    truepivsB(_,k) = truepivsB_k;
  }
  List ret;
  ret["truepivsA"] = truepivsA;
  ret["truepivsB"] = truepivsB;
  return ret;
}