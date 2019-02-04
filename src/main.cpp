#include <unistd.h>
#include <iostream>

#include "mmseqs.h"
#include "mafft.h"
#include "format.h"
#include "nj.h"
#include "transitivity.h"

using namespace std;

int main(int argc, char* argv[]){
  
  /*/Check the mmseqs&mafft command/*/
  auto status1 = system("which mmseqs &> /dev/null");
  auto status2 = system("which mafft  &> /dev/null");
  if(WEXITSTATUS(status1) != 0){
    return -1;
  }
  if(WEXITSTATUS(status2) != 0){
    return -1;
  }

  /*/Getopt/*/
  //int silence = 0;              // -s
  //int ep_num  = 0;              // -e
  //int seed    = 0;              // -r
  //int label   = 0;              // -l
  string threads = "1";           // -t
  string sensitivity = "7.5";     // -m
  string deletion_method = "pgd"; // -g [default: pgd=pairwise gap deletion]
  double wp     = 0.5;            // -w
  //string bs_method = "";        // -b [fbp (Felsenstein's bootstrap proportion) or tbe (transfer bootstrap expectation)]

  opterr = 0; // default error messages -> OFF
  int opt;
  while ((opt = getopt(argc, argv, "g:w:")) != -1){
    if(opt == 'w'){ // OK! (./gs -w 0.5 IN.fst)
      wp = atof(optarg);
    }
    else if(opt == 'g'){ // OK! (./gs -w 0.5 IN.fst)
      deletion_method = optarg;
    }
  }

  /*/Input file name/*/
  string input = "";
  if(optind < argc){ // OK! (./gs -w 0.5 IN.fst)
    input = argv[optind];
  }
  else{ // NG! (./gs -w 0.5)
    return -1;
  }
  
  /*/Variables/*/
  double* Wp;      // Distance matrix
  double* Wm;      // Distance matrix
  double* Wpm;     // Distance matrix
  double* M;       // Sequence similarity matrix
  int size;        // Row size of W (W is a synmetry matrix)
  int* njpm;       // Result of NJ method
  string newickp;  // NJ tree (without EP values)
  string newickm;  // NJ tree (without EP values)
  string newickpm; // NJ tree (without EP values)

  /*/File I/O/*/
  regex re(R"(\.[^\.]+$)"); // Extention of the input file
  auto original_fasta = string(input);
  auto annotation_txt = regex_replace(original_fasta, re, "_annotation.txt");
  auto simple_fasta   = regex_replace(original_fasta, re, "_simple.fst");
  auto mmseqs_result  = regex_replace(original_fasta, re, "_mmseqs.txt");
  auto mafft_result   = regex_replace(original_fasta, re, "_mafft.txt");
  
  ifstream ifs1(original_fasta); // Fasta file (original)
  ofstream ofs1(annotation_txt); // Annotation file
  ofstream ofs2(simple_fasta);   // Fasta file (simple)
  if(ifs1.fail()){
    return -1;
  }
  if(ofs1.fail()){
    return -1;
  }
  if(ofs2.fail()){
    return -1;
  }

  /*/Parsing fasta file/*/
  auto file_condition = readFASTA(ifs1, ofs1, ofs2, size);
  if(file_condition == 1){
    return -1;
  }
  else if(file_condition == 2){
    return -1;
  }
  else if(file_condition == 3){
    return -1;
  }
  else if(file_condition == 4){
    return -1;
  }  
  
  /*/Executing MMSeqs/*/
  mmseqs(simple_fasta, mmseqs_result, threads, sensitivity);
  ifstream ifs2(mmseqs_result); // MMseqs result file
  if(ifs2.fail()){
    return -1;
  }

  /*/Executing MAFFT/*/
  mafft(simple_fasta, mafft_result);
  ifstream ifs3(mafft_result); // MAFFT result file
  if(ifs3.fail()){
    return -1;
  }

  /*/Reading Data/*/
  auto same_sequence = 
  PSA2mat(ifs2, Wp, M, size, deletion_method);
  MSA2mat(ifs3, Wm, size, deletion_method);
  if(same_sequence>0){}

  //auto transitivity_score = transitivity(M, size);

  /*/Mixed Inference/*/
  mix(Wp, Wm, Wpm, wp, size);

  /*/NJ method/*/
  NJ(Wpm, njpm, size);

  /*/Generating NJ tree Newick/*/
  sc2nwk(njpm, newickpm, size);

  /*/GS tree WITHOUT EP values ->STDOUT/*/
  cout << newickpm << endl;

  return 0;
}
