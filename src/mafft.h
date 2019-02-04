/******************************************\
| Graph Splitting Method v2.2 (2018/11/16) |
|                                          |
| Copyright (c) 2015-2018 Motomu Matsui    |
|     Distributed under the GNU GPL        |
|                                          |
|     Matsui M and Iwasaki W (2018)        |
|     Systematic Biology, xx:xx-xx.        |
|                                          |
|     http://gs.bs.s.u-tokyo.ac.jp/        |
\******************************************/

#ifndef MAFFT_H
#define MAFFT_H

#include <fstream>
#include <regex>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unordered_map>

using namespace std;

void mafft(string const& fst, string const& output);

#endif /*MAFFT_H*/
