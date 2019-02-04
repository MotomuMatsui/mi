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

#ifndef DISTANCE_H
#define DISTANCE_H

#include <iostream>
#include <string>
#include <cmath>

using namespace std;

double sum_score(string const& A, string const& B, string const& opt);
double scoredist(string const& A, string const& B, string const& opt);
double scoredist2(
		  double const& xx,  double const& yy,  double const& xy,  double const& yx, 
		  double const& lxx, double const& lyy, double const& lxy, double const& lyx);

#endif /*DISTANCE_H*/
