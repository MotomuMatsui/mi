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

#include <fstream>
#include <regex>
#include <sstream>
#include <stdio.h>
#include <string>
#include <sys/stat.h>
#include <unordered_map>

#include "mafft.h"

using namespace std;

void mafft(string const& fst, string const& output){
  //mmseqs command
  const string mafft = "mafft";

  //Commands
  auto cmd  = mafft+" --quiet --anysymbol "+fst+" > "+output;
  auto info = system(cmd.c_str()); if(info>0){}
}
