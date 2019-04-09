#include <random>   
#include <regex>
#include <unistd.h>
#include <iostream>

#include "ep.h"
#include "nj.h"
#include "mmseqs.h"
#include "mafft.h"
#include "messages.h"
#include "format.h"
#include "transitivity.h"

using namespace std;

int main(int argc, char* argv[]){
  
  /*/Check the mmseqs&mafft command/*/
  auto status1 = system("which mmseqs &> /dev/null");
  auto status2 = system("which mafft  &> /dev/null");
  if(WEXITSTATUS(status1) != 0){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "mmseqs not found!\nCheck our web page (https://github.com/MotomuMatsui/mi) for more information" << endl;

    return -1;
  }
  if(WEXITSTATUS(status2) != 0){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "mafft not found!\nCheck our web page (https://github.com/MotomuMatsui/mi) for more information" << endl;

    return -1;
  }

  /*/Getopt/*/
  int silence            = 0;     // -s
  int ep_num             = 0;     // -e
  int seed               = 0;     // -r
  int label              = 0;     // -l
  string threads         = "1";   // -t
  string sensitivity     = "7.5"; // -m
  string deletion_method = "cgd"; // -g [default: cgd=complete gap deletion]
  double wp              = 0.8;   // -w
  string bs_method       = "";    // -b [fbp (Felsenstein's bootstrap proportion) or tbe (transfer bootstrap expectation)]

  opterr = 0; // default error messages -> OFF
  int opt;
  regex renum(R"(^[\d\.]+$)"); // -e/-r/-t/-m option requires an integer/flout number
  while ((opt = getopt(argc, argv, "shlve:r:t:m:b:g:w:")) != -1){
    if(opt == 'e'){ // OK! (./gs -e 100 IN.fst)
      if(regex_match(optarg, renum)){
        ep_num = atoi(optarg);
      }
      else{ // NG! (./gs -e hundred IN.fst)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
    else if(opt == 'r'){ // OK! (./gs -r 12345 IN.fst)
      if(regex_match(optarg, renum)){
        seed = atoi(optarg);
      }
      else{ // NG! (./gs -r one_two_three IN.fst)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -r requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
    else if(opt == 't'){ // OK! (./gs -t 4 IN.fst)
      if(regex_match(optarg, renum)){
        auto th = atoi(optarg);
        if(th >0){
          threads = string(optarg);
        }
        else{ // NG! (./gs -t four IN.fst)
          /*PRINT*/ print_banner();
          /*PRINT*/ cerr << "Option -t requires an integer argument (>=1).\n" << endl;
          /*PRINT*/ print_usage(argv[0]);
          return -1;
        }
      }
      else{ // NG! (./gs -t four IN.fst)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -t requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
    else if(opt == 'm'){ // OK! (./gs -m 7.5 IN.fst)
      if(regex_match(optarg, renum)){
        auto sen = atof(optarg);
        if(1<=sen && sen<=7.5){
          sensitivity = string(optarg);
        }
        else{ // NG! (./gs -m 10 IN.fst)
          /*PRINT*/ print_banner();
          /*PRINT*/ cerr << "Option -m requires a double number argument [1, 7.5].\n" << endl;
          /*PRINT*/ print_usage(argv[0]);
          return -1;
        }
      }
      else{ // NG! (./gs -m seven IN.fst)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -m requires a double number argument [1, 7.5].\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
    else if(opt == 'h'){ // HELP message (./gs -h)
      /*PRINT*/ print_banner();
      /*PRINT*/ print_usage(argv[0]);
      return 0;
    }
    else if(opt == 'v'){ // Version (./gs -v)
      /*PRINT*/ print_banner();
      return 0;
    }
    else if(opt == 's'){ // SILENT mode (./gs -s -e 100 IN.fst)
      silence = 1;
    }
    else if(opt == 'l'){ // LABEL mode (./gs -l -e 100 IN.fst)
      label = 1;
    }
    else if(opt == 'b'){ // Statistical method to assess the robustness of inffered branches (./gs -e 100 -b tbe IN.fst)
      bs_method = optarg;
      if(bs_method != "fbs" && bs_method != "tbe"){
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -b requires a string that equals 'fbs' or 'tbe'.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
    else if(opt == 'w'){ // OK! (./gs -w 0.5 IN.fst)
      wp = atof(optarg);
    }
    else if(opt == 'g'){ // OK! (./gs -w 0.5 IN.fst)
      deletion_method = optarg;
    }
    else if (opt == '?'){
      if(optopt == 'e'){ // NG! (./gs IN.fst -e)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -e requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
      else if(optopt == 'r'){ // NG! (./gs IN.fst -r)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -r requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
      else if(optopt == 't'){ // NG! (./gs IN.fst -t)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -t requires an integer argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
      else if(optopt == 'm'){ // NG! (./gs IN.fst -m)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -m requires a flout argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
      else if(optopt == 'b'){ // NG! (./gs IN.fst -b)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << "Option -b requires a string argument.\n" << endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
      else{ // NG! (./gs -Z)
        /*PRINT*/ print_banner();
        /*PRINT*/ cerr << argv[0] << ": invalid option\n" <<  endl;
        /*PRINT*/ print_usage(argv[0]);
        return -1;
      }
    }
  }

  /*/Input file name/*/
  string input = "";
  if(optind < argc){ // OK! (./gs -w 0.5 IN.fst)
    input = argv[optind];
  }
  else{ // NG! (./gs -w 0.5)
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << argv[0] << " requires an input file (fasta format).\n" << endl;
    /*PRINT*/ print_usage(argv[0]);
    return -1;
  }
  
  /*/Variables/*/
  double* Wp;      // Distance matrix
  double* Wm;      // Distance matrix
  double* Wpm;     // Distance matrix
  double* M;       // Sequence similarity matrix
  int size;        // Row size of W (W is a synmetry matrix)
  int* njpm;       // Result of NJ method
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
    /*PRINT*/ cerr << "\nCannot access " << original_fasta << "!" << endl;
    return -1;
  }
  if(ofs1.fail()){
    /*PRINT*/ cerr << "\nCannot create " << annotation_txt << "!" << endl;
    return -1;
  }
  if(ofs2.fail()){
    /*PRINT*/ cerr << "\nCannot create " << simple_fasta << "!" << endl;
    return -1;
  }

  /*/Parsing fasta file/*/
  auto file_condition = readFASTA(ifs1, ofs1, ofs2, size);
  ofs1.close();
  ofs2.close();

  if(file_condition == 1){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "Input file had an empty entry.\n" << endl;
    /*PRINT*/ print_usage(argv[0]);    
    return -1;
  }
  else if(file_condition == 2){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "Input file had an empty entry (1st entry lacked '>').\n" << endl;
    /*PRINT*/ print_usage(argv[0]);
    return -1;
  }
  else if(file_condition == 3){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "Input file had an empty entry (last entry lacked '>').\n" << endl;
    /*PRINT*/ print_usage(argv[0]);    
    return -1;
  }
  else if(file_condition == 4){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "Input file have to contain more than two sequences.\n" << endl;
    /*PRINT*/ print_usage(argv[0]);    
    return -1;
  }  

  
  /*/ Parameters /*/  
  if(!silence){
    /*PRINT*/ print_banner();
    /*PRINT*/ cerr << "Settings:" << endl;
    /*PRINT*/ cerr << "-Input" << endl;
    /*PRINT*/ cerr << "  file = " << input << endl;
    /*PRINT*/ cerr << "  # of sequences = " << size << endl << endl;

    /*PRINT*/ cerr << "-MMseqs" << endl;
    /*PRINT*/ cerr << "  sentitivity = " << sensitivity << endl;
    /*PRINT*/ cerr << "  # of threads = " << threads << endl << endl;

    /*PRINT*/ cerr << "-EP method" << endl;
    if(ep_num>0){
      if(bs_method == "fbs"){
        /*PRINT*/ if(!silence) cerr << "  method = Felsenstein's bootstrap proportion" << endl;
      }
      else{
        /*PRINT*/ if(!silence) cerr << "  method = Transfer bootstrap expectation (F. Lemoine, et al., Nature, 2018)" << endl;
      }
    }

    if(seed>0){
      /*PRINT*/ cerr << "  random seed = " << seed << endl;
    }
    else{
      /*PRINT*/ cerr << "  random seed = " << "a random number (default)" << endl;
    }
    /*PRINT*/ cerr << "  # of iterations = " << ep_num << endl << endl;
    /*PRINT*/ cerr << "Progress:" << endl;
  }
  
  /*/Executing MMSeqs/*/
  /*PRINT*/ if(!silence) cerr << "-MMseqs\n  " << size << "x" << size << " pairwise alignment\n" << "  searching...\r" << flush;
  mmseqs(simple_fasta, mmseqs_result, threads, sensitivity);
  ifstream ifs2(mmseqs_result); // MMseqs result file
  if(ifs2.fail()){
    /*PRINT*/ cerr << "  Cannot access " << mmseqs_result << "!" << endl;
    return -1;
  }
  else{
    /*PRINT*/ if(!silence) cerr << "  done.        " << endl << endl;
  }

  /*/Executing MAFFT/*/
  /*PRINT*/ if(!silence) cerr << "-Mafft\n  " << size << "x" << size << " multiple alignment\n" << "  searching...\r" << flush;
  mafft(simple_fasta, mafft_result);
  ifstream ifs3(mafft_result); // MAFFT result file
  if(ifs3.fail()){
    /*PRINT*/ cerr << "  Cannot access " << mafft_result << "!" << endl;
    return -1;
  }
  else{
    /*PRINT*/ if(!silence) cerr << "  done.        " << endl << endl;
  }

  /*/Reading Data/*/
  auto same_sequence = 
  PSA2mat(ifs2, Wp, M, size, deletion_method);
  MSA2mat(ifs3, Wm, size, deletion_method);  
  /*PRINT*/ if(!silence) if(same_sequence>0) cerr << "  <WARNING> This dataset has " << same_sequence << " duplicated sequence pair(s)" << endl << endl;

  auto transitivity_score = transitivity(M, size);
  /*PRINT*/ if(!silence) cerr << "  transitivity = "<< transitivity_score << endl << endl;

  /*/Mixed Inference/*/
  mix(Wp, Wm, Wpm, wp, size);

  /*/NJ method/*/
  /*PRINT*/ if(!silence) cerr << "-NJ method\n" << "  executing...\r" << flush;
  NJ(Wpm, njpm, size);
  /*PRINT*/ if(!silence) cerr << "  done.         " << endl << endl;

  /*/Generating NJ tree Newick/*/
  sc2nwk(njpm, newickpm, size);

  /*/ EP method /*/
  if(ep_num>0){
    /*PRINT*/ if(!silence) cerr << "-EP method" << endl;

    unordered_map<string, double> ep;
    string newick_EP; // GS+EP tree

    if(bs_method == "fbs"){
      // Random number generator (Uniform distribution->Mersenne Twister)
      function<double()> R;
      uniform_real_distribution<double> urd(0,1);    // uniform distributed random number
      
      if(seed>0){
        mt19937 mt(static_cast<unsigned int>(seed)); // mersenne twister
        R = bind(urd, ref(mt));                      // random number generator    
      }
      else{
        random_device rd;                            // random seed
        mt19937 mt(rd());                            // mersenne twister
        R = bind(urd, ref(mt));                      // random number generator        
      }    

      for(int n=1; n<=ep_num; n++){
        /*PRINT*/ if(!silence) cerr << "  " << n << "/" << ep_num << " iterations" << "\r"<< flush;

        EP_fbs(Wpm, ep, R, size);
	// W: INPUT (sequence similarity matrix)
	// ep: OUTPUT (result of Edge Perturbation method)
	// R: random number generator
      }
    }
    else{
      int* list_ori; 
      sc2list(njpm, list_ori, size);
      // gs: INPUT (result of stepwise spectral clustering)
      // list: OUTPUT (NJ tree [leaves])      

      // Random number generator (Uniform distribution->Mersenne Twister)
      function<double()> R;
      uniform_real_distribution<double> urd(0,1);    // uniform distributed random number
      
      if(seed>0){
        mt19937 mt(static_cast<unsigned int>(seed)); // mersenne twister
        R = bind(urd, ref(mt));                      // random number generator    
      }
      else{
        random_device rd;                            // random seed
        mt19937 mt(rd());                            // mersenne twister
        R = bind(urd, ref(mt));                      // random number generator        
      }    

      for(int n=1; n<=ep_num; n++){
        /*PRINT*/ if(!silence) cerr << "  " << n << "/" << ep_num << " iterations" << "\r"<< flush;

        EP_tbe(Wpm, list_ori, ep, R, size);
	// W: INPUT (sequence similarity matrix)
	// ep: OUTPUT (result of Edge Perturbation method)
	// R: random number generator
      }

      delete[] list_ori;
    }

    /*PRINT*/ if(!silence) cerr << "\n  done." << endl << endl;
    /*PRINT*/ if(!silence) cerr << "------------------------------------------\n" << endl;

    addEP(newickpm, newick_EP, ep, ep_num, size);
    // newick: INPUT (GS tree [newick format])
    // newick_EP: OUTPUT (GS+EP tree [newick format])
    // ep: INPUT (result of Edge Perturbation method)
    // ep_num: INPUT (# of Edge Perturbation method)

    /*/ GS tree WITH EP values ->STDOUT /*/
    if(label==1){
      string newick_ann;
      addLABEL(newick_EP, newick_ann, annotation_txt, size);

      cout << newick_ann << endl;
    }
    else{
      cout << newick_EP << endl;
    }
  }
  else{ // skip the EP method
    /*PRINT*/ if(!silence) cerr << "------------------------------------------\n" << endl;

    /*/ GS tree WITHOUT EP values ->STDOUT /*/
    if(label==1){
      string newick_ann;
      addLABEL(newickpm, newick_ann, annotation_txt, size);

      cout << newick_ann << endl;
    }
    else{
      cout << newickpm << endl;
    }
  }

  delete[] Wpm;
  delete[] njpm;

  return 0;
}
