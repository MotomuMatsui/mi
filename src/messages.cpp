#include <iostream>
#include <string>

#include "messages.h"

using namespace std;

void print_banner(){

  string banner = 
    "------------------------------------------\n"
    " Mixed Inference Method v1.0 (2019/04/10) \n"
    "                                          \n"
    "   Copyright (c) 2019 Motomu Matsui       \n"
    "   ???, xx:xx-xxx, 2019                   \n"
    "                                          \n"
    "   http://gs.bs.s.u-tokyo.ac.jp/          \n"
    "------------------------------------------\n";

  cerr << banner << endl;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////

void print_usage(char*& program){
  cerr << "Usage: " << program << " [-e INTEGER(>=0)] [-b STRING(tbe/fbs)] [-r INTEGER(>0)] [-t INTEGAR(>0)] [-m FLOAT(1-7.5)] [-s] [-l] [-h] [-v] input > output" << endl;
  cerr << "-e " << "the number of replicates for EP method. Default: 0" << endl;
  cerr << "-b " << "the bootstrap method, tbe or fbs. Default: tbe" << endl;
  cerr << "-r " << "the random seed number for EP method. Default: random number" << endl;
  cerr << "-t " << "the number of threads for MMseqs. Default: 1" << endl;
  cerr << "-m " << "sensitivity for MMseqs. Default: 7.5" << endl;
  cerr << "-s " << "silent mode: do not report progress. Default: Off" << endl;
  cerr << "-l " << "newick format with actual names. Default: Off" << endl;
  cerr << "-h " << "show help messages. Default: Off" << endl;
  cerr << "-v " << "show the version. Default: Off" << endl;
}

//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
