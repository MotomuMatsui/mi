#include "ep_function.h"

using namespace std;

//Generalized Extreme Value function (inverse function)                                                                        
double gev(double const& x, double const& mu){
  double theta = mu*(1-mu)/3;
  double gamma = exp(-3*mu)-1;

  if(gamma == 0){ //Gummbel distribution                                                                                        
    return mu - theta*log(-log(x));
  }
  else{
    return mu +( pow(-log(x),-gamma)-1 )*theta/gamma;
  }
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
