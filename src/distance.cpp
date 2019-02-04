/********************************************\
| Edge Perturbation Method v2.0 (2018/11/16) |
|                                            |
|  Copyright (c) 2015-2018 Motomu Matsui     |
|      Distributed under the GNU GPL         |
|                                            |
|      Matsui M and Iwasaki W (2018)         |
|      Systematic Biology, xx:xx-xx.         |
|                                            |
|      http://gs.bs.s.u-tokyo.ac.jp/         |
\********************************************/

#include "distance.h"
#include "matrix.h"
#include <iostream>
#include <string>
#include <cmath>
#include <unistd.h>

using namespace std;

double scoredist(string const& A, string const& B, string const& opt){
  
  /*Parameters*/
  int sc   = 0;
  int scA  = 0;
  int scB  = 0;
  int len  = A.size();
  int lenA = 0;
  int lenB = 0;

  if(opt == "pgd"){
    for(int p = 0; p < len; p++){
      // gap; sc -= 0.6;      
      
      int pa = s2p[(unsigned char)(A[p])];
      int pb = s2p[(unsigned char)(B[p])];
      
      sc  += BLOSUM62p[pa][pb];
      scA += BLOSUM62p[pa][pa];
      scB += BLOSUM62p[pb][pb];

      lenA += PairwiseGapDeletion[pa][pa];
      lenB += PairwiseGapDeletion[pb][pb];
    }
  }
  else{
    for(int p = 0; p < len; p++){
      // gap; sc -= 0.6;      
      
      int pa = s2p[(unsigned char)(A[p])];
      int pb = s2p[(unsigned char)(B[p])];
      
      sc  += BLOSUM62c[pa][pb];
      scA += BLOSUM62c[pa][pa];
      scB += BLOSUM62c[pb][pb];

      lenA += CompleteGapDeletion[pa][pa];
      lenB += CompleteGapDeletion[pb][pb];
    }
  }

  double scMAX = (scA + scB)/2;
  double scR   = -0.5209 * (lenA+lenB)/2;
  double od    = (sc - scR)/(scMAX - scR);
  
  if (od < 0.05) od = 0.01;  /* Limit to 300 PAM;  len==0 if no overlap */
  if (od > 1.0)  od = 1.0;
  
  //double cd = -log(od) * 100.0 * 1.337;
  double cd = -log(od)*1.2110;

  return cd;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

double sum_score(string const& A, string const& B, string const& opt){
  
  /*Parameters*/
  int len = A.size();
  int sc  = 0;

  if(opt == "pgd"){
    for(int p = 0; p < len; p++){
      // gap; sc -= 0.6;
      
      int pa = s2p[(unsigned char)(A[p])];
      int pb = s2p[(unsigned char)(B[p])];
      
      sc += BLOSUM62p[pa][pb];
    }
  }
  else{
    for(int p = 0; p < len; p++){
      // gap; sc -= 0.6;
      
      int pa = s2p[(unsigned char)(A[p])];
      int pb = s2p[(unsigned char)(B[p])];
      
      sc += BLOSUM62c[pa][pb];
    }
  }

  return sc;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

double scoredist2(
		  double const& xx,  double const& yy,  double const& xy,  double const& yx, 
		  double const& lxx, double const& lyy, double const& lxy, double const& lyx
){
  double scMEAN  = (xy + yx)/2;
  double scMAX   = (xx + yy)/2;
  double scRmean = -0.5209 * (lxy+lyx)/2;
  double scRmax  = -0.5209 * (lxx+lyy)/2;
  double od      = (scMEAN - scRmean)/(scMAX - scRmax);
  
  if (od < 0.01) od = 0.01;  /* Limit to 300 PAM;  len==0 if no overlap */
  if (od > 1.0)  od = 1.0; 
    
  double cd = -log(od)*1.2110; //original:1.337
  
  return cd;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

