/******************************************\
| Graph Splitting Method v2.3 (2018/11/16) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#ifndef FORAMT_H
#define FORMAT_H

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <regex>

using namespace std;

int readFASTA(ifstream& ifs, ofstream& ofs1, ofstream& ofs2, int& row);
int PSA2mat(ifstream& ifs, double* (&W), double* (&M), int const& size, string const& deletion_method);
int MSA2mat(ifstream& ifs, double* (&W), int const& size, string const& deletion_method);
void sc2nwk(int* const& W, string& newick, int const& size);
void addEP(string const& newick, string& newick_EP, unordered_map<string, double>& ep, int const& ep_num, int const& size);
void addLABEL(string const& newick, string& newick_ann, string const& annotation_txt, int const& size);
void sc2list(int* const& gs, int* (&list), int const& size);
void mix(double* const& Wp, double* const& Wm, double* (&Wpm), double const& wp, int const& size);

#endif /*FORMAT_H*/
