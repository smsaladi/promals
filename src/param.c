#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include "param.h"

double ave_grp_thr = 1;
double minProb = 0.01;
string outFile = "";
int probconsBLOSUM = 1;  // probcons blosum62
int useLocal = 3;
float weightG = 1;
int ss = 3;
int solv = 1;
int unaligned = 0;
char parameter_file[200];
char parameter_file1[200];
char parameter_file2[200];
int relax_number = 2;
int reverse_align_order = 0;
double id_thr = 0.6;

char program_dir[500];

int psipred_env_number;
char psipred_dir[500];
char runpsipred_command[500];
char runpsipred1_command[500];

// for database homolog
char blast_dir[500];
char uniref90_file[500];
char blastpgp_cmd[500];
char blastpgp_command[500];

int use_ss_freq = 2;

char psipred_parameter_file[200];
int use_single_sequence = 0;

// float ss_w = 0.3;
// float score_w = 1;
float ss_w = 0.2;
float score_w = 0.8;
float score_shift = 0;
int adjust_weight = 0;

int relax_count = 0;

// blastpgp options
int iteration_number = 3;
double evalue = 0.001;
char filter_query[3];

// purge blast output option
int max_num_sequences = 300;
double low_id_thr = 0.2;

// output alignment
int blocksize = 70;

// maximum number of pre-aligned groups
int max_group_number = 60;

// clean blastpgp and psipred intermediate results
int clean_blast_after = 0;
int clean_blast_before = 0;

// PSIBLAST alignment of structures: identity cutoff
int struct_id_cutoff = 30;
int below_id_cutoff = 101;

// before or after relax: combine structure information
int before_relax_combine = 1;  // default: combine after two relaxation rounds

// PSIBLAST alignment of structures: identity cutoff
float struct_weight = 1.5;  // derived from structural homologs
float sequence_weight =
    1;  // derived from profile-profile alignment with secondary structure
float user_constraint_weight = 1.5;  // derived from user-defined constraints
float pdb_weight = 1.5;  // derived from structure alignment of input structures
double minDaliZ = 1.0;

// use different structural alignments
int use_dali = 0;
int use_fast = 1;
int use_tmalign = 1;

// realign psiblast alignment
int realign_psiblast = 0;

// constraint file
char constraint_file[500];  // for structural alignments
char user_constraint[500];  // for user-defined constraints

//
double weight_end_penalty = 0.75;

//
char mafft[500];

// filter similar
int filter_similar = 1;
char cdhit[500];
float cdhit_c_option = 0.95;
int exclude_similar = 0;

// updated database
int use_updated_database = 1;

//
int max_struct_number = 1;

void getParameter(int argc, char **argv, int prog) {
  int i, j;
  string argStr;

  if (argc <= 1) {
    printHelp(prog);
    exit(0);
  }

  if (argv[1][0] == '-') {
    printHelp(prog);
    exit(0);
  }

  // cout << argc << endl;
  parameter_file[0] = '\0';
  psipred_dir[0] = '.';
  psipred_dir[1] = '\0';
  blast_dir[0] = '.';
  blast_dir[1] = '\0';
  blastpgp_cmd[0] = '\0';
  uniref90_file[0] = '\0';
  psipred_parameter_file[0] = '\0';
  psipred_env_number = 3;

  strcpy(filter_query, "T\0");

  constraint_file[0] = '\0';
  user_constraint[0] = '\0';

  strcpy(blastpgp_command, "blastpgp");
  strcpy(runpsipred_command, "runpsipredplus");
  strcpy(runpsipred1_command, "runpsipredplus");

  if (getenv("PROMALS_DIR")) {
    strcpy(program_dir, getenv("PROMALS_DIR"));
    strcat(program_dir, "/");
  } else {
    cout << "Please set PROMALS_DIR. Otherwise, set to `.`" << endl;
    strcpy(program_dir, "./");
  }

  strcpy(uniref90_file, program_dir);
  strcat(uniref90_file, "/db/uniref90/uniref90_filt");

  strcpy(psipred_parameter_file, program_dir);
  strcat(psipred_parameter_file, "/db/param/dataset_0.20_0.50_0.60_abcd.mat");

  strcpy(mafft, "mafft");
  strcpy(cdhit, "cd-hit");

  strcpy(blast_dir, argv[1]);
  strcat(blast_dir, "_blast");

  strcpy(parameter_file,
         "hmm_parameters/dataset_0.20_0.40_0.60_abcd.dali0.solv1_ss1.mat");

  for (i = 1; i < argc; i++) {
    argStr = argv[i];
    // if(argStr == "-a") { ave_grp_thr = atof(argv[i+1]); continue; }
    if (argStr == "-minprob") {
      minProb = atof(argv[i + 1]);
      continue;
    }
    if (argStr == "-outfile") {
      outFile = argv[i + 1];
      continue;
    }
    if (argStr == "-pb") {
      probconsBLOSUM = atoi(argv[i + 1]);
      continue;
    }
    if (argStr == "-local") {
      useLocal = atoi(argv[i + 1]);
    }
    if (argStr == "-wg") {
      weightG = atof(argv[i + 1]);
      assert(weightG >= 0);
      assert(weightG <= 1);
    }
    if (argStr == "-ss") {
      ss = atoi(argv[i + 1]);
      assert(ss > 0);
      assert(ss <= 3);
      assert(ss != 2);
    }
    if (argStr == "-solv") {
      solv = atoi(argv[i + 1]);
      assert(solv > 0);
      assert(solv <= 3);
    }
    if (argStr == "-unaligned") {
      unaligned = atoi(argv[i + 1]);
      assert(unaligned <= 1);
    }

    if (argStr == "-param") {
      strcpy(parameter_file, argv[i + 1]);
    }
    if (argStr == "-param1") {
      strcpy(parameter_file1, argv[i + 1]);
    }
    if (argStr == "-param2") {
      strcpy(parameter_file2, argv[i + 1]);
    }
    if (argStr == "-r") {
      relax_number = atoi(argv[i + 1]);
    }
    if (argStr == "-reverse") {
      reverse_align_order = atoi(argv[i + 1]);
    }
    if (argStr == "-id_thr") {
      id_thr = atof(argv[i + 1]);
      // cout << id_thr << endl;
    }

    if (argStr == "-psipred_dir") {
      strcpy(psipred_dir, argv[i + 1]);
    }
    if (argStr == "-env_number") {
      psipred_env_number = atoi(argv[i + 1]);
    }
    if (argStr == "-psipred_param") {
      strcpy(psipred_parameter_file, argv[i + 1]);
    }
    if (argStr == "-use_single_sequence") {
      use_single_sequence = atoi(argv[i + 1]);
    }
    if (argStr == "-ss_weight") {
      ss_w = atof(argv[i + 1]);
    }
    if (argStr == "-score_weight") {
      score_w = atof(argv[i + 1]);
    }
    if (argStr == "-score_shift") {
      score_shift = atof(argv[i + 1]);
    }
    if (argStr == "-blast_dir") {
      strcpy(blast_dir, argv[i + 1]);
    }
    if (argStr == "-use_ss_freq") {
      use_ss_freq = atoi(argv[i + 1]);
    }
    if (argStr == "-adjust_weight") {
      adjust_weight = atoi(argv[i + 1]);
    }
    if (argStr == "-relax_count") {
      relax_count = atoi(argv[i + 1]);
    }
    if (argStr == "-iter_number") {
      iteration_number = atoi(argv[i + 1]);
      if (iteration_number <= 0) {
        cout << "Error: iteration number must be a positive integer" << endl;
        exit(0);
      }
      iteration_number++;  // to get a check point file, iteration number must
                           // be above 1
      // if(iteration_number <=0) iteration_number = 1;
      if (iteration_number > 50) iteration_number = 50;
    }
    if (argStr == "-evalue") {
      evalue = atof(argv[i + 1]);
      if (evalue <= 0) {
        cout << "Error: evalue must be a positve value" << endl;
        exit(0);
      }
    }
    if (argStr == "-filter_query") {
      strncpy(filter_query, argv[i + 1], 1);
      if (filter_query[0] == 'f') filter_query[0] = 'F';
      if (filter_query[0] != 'F') filter_query[0] = 'T';
    }
    if (argStr == "-max_num_sequences") {
      max_num_sequences = atoi(argv[i + 1]);
      if (max_num_sequences <= 0) {
        cout << "Error: number of selected homologs must be positive" << endl;
        exit(0);
      }
    }
    if (argStr == "-low_id_thr") {
      low_id_thr = atof(argv[i + 1]);
      if (low_id_thr > 1) {
        cout << "Error: identity cutoff for excluding divergent homologs must "
                "be less than 1"
             << endl;
        exit(0);
      }
    }
    if (argStr == "-blocksize") {
      blocksize = atoi(argv[i + 1]);
      if (blocksize < 10) {
        blocksize = 70;
      }
    }
    if (argStr == "-max_group_number") {
      max_group_number = atoi(argv[i + 1]);
      if (max_group_number < 1) {
        max_group_number = 1;
      }
    }
    if (argStr == "-clean_blast_after") {
      clean_blast_after = atoi(argv[i + 1]);
    }
    if (argStr == "-clean_blast_before") {
      clean_blast_before = atoi(argv[i + 1]);
    }
    if (argStr == "-struct_id_cutoff") {
      struct_id_cutoff = atoi(argv[i + 1]);
      if (struct_id_cutoff < 1) {
        struct_id_cutoff = 1;
      }
    }
    if (argStr == "-below_id_cutoff") {
      below_id_cutoff = atoi(argv[i + 1]);
      if (below_id_cutoff < 1) {
        below_id_cutoff = 1;
      }
    }
    if (argStr == "-before_relax_combine") {
      before_relax_combine = atoi(argv[i + 1]);
      if (before_relax_combine < 0) {
        before_relax_combine = 0;
      }
      if (before_relax_combine > 2) {
        before_relax_combine = 2;
      }
    }
    if (argStr == "-struct_weight") {
      struct_weight = atof(argv[i + 1]);
      if (struct_weight < 0) {
        struct_weight = 0;
      }
    }
    if (argStr == "-sequence_weight") {
      sequence_weight = atof(argv[i + 1]);
      if (sequence_weight < 0) {
        sequence_weight = 0;
      }
    }
    if (argStr == "-user_constraint_weight") {
      user_constraint_weight = atof(argv[i + 1]);
      if (user_constraint_weight < 0) {
        user_constraint_weight = 0;
      }
    }
    if (argStr == "-pdb_weight") {
      pdb_weight = atof(argv[i + 1]);
      if (pdb_weight < 0) {
        pdb_weight = 0;
      }
    }
    if (argStr == "-minDaliZ") {
      minDaliZ = atof(argv[i + 1]);
      if (minDaliZ < 0) {
        minDaliZ = 0;
      }
    }
    if (argStr == "-dali") {
      use_dali = atoi(argv[i + 1]);
      // if(dali < 0) { dali = 0; }
      if (use_dali) {
        char command1[500];
        sprintf(command1, "ls %s/bin/DaliLite", program_dir);
        int a1 = system(command1);
        if (a1 != 0) {
          cout << "Error: DaliLite program is not available at " << program_dir
               << endl;
          exit(1);
        }
      }
    }
    if (argStr == "-fast") {
      use_fast = atoi(argv[i + 1]);
      // if(minDaliZ < 0) { minDaliZ = 0; }
      if (use_fast) {
        char command1[500];
        sprintf(command1, "ls %s/bin/fast", program_dir);
        int a1 = system(command1);
        if (a1 != 0) {
          cout << "Error: fast program is not available at " << program_dir
               << endl;
          exit(1);
        }
      }
    }
    if (argStr == "-tmalign") {
      use_tmalign = atoi(argv[i + 1]);
      // if(minDaliZ < 0) { minDaliZ = 0; }
      if (use_tmalign) {
        char command1[500];
        sprintf(command1, "ls %s/bin/tmalign", program_dir);
        int a1 = system(command1);
        if (a1 != 0) {
          cout << "Error: TMalign executable file is not available at "
               << program_dir << endl;
          exit(1);
        }
      }
    }
    if (argStr == "-realign_psiblast") {
      realign_psiblast = atoi(argv[i + 1]);
      // if(minDaliZ < 0) { minDaliZ = 0; }
    }
    if (argStr == "-constraint") {
      strcpy(constraint_file, argv[i + 1]);
      // if(minDaliZ < 0) { minDaliZ = 0; }
    }
    if (argStr == "-user_constraint") {
      strcpy(user_constraint, argv[i + 1]);
      // if(minDaliZ < 0) { minDaliZ = 0; }
    }
    if (argStr == "-weight_end_penalty") {
      weight_end_penalty = atof(argv[i + 1]);
      if (weight_end_penalty < 0) {
        weight_end_penalty = 1.0;
      }
    }
    if (argStr == "-filter_similar") {
      filter_similar = atoi(argv[i + 1]);
    }
    if (argStr == "-exclude_similar") {
      exclude_similar = atoi(argv[i + 1]);
    }
    if (argStr == "-use_updated_database") {
      use_updated_database = atoi(argv[i + 1]);
    }
    if (argStr == "-max_struct_number") {
      max_struct_number = atoi(argv[i + 1]);
    }
  }
  // cout << "Here" << endl;

  // check if blast_dir is a regular file
  struct stat buf;

  if (lstat(blast_dir, &buf) < 0) {
    cout << "blast_dir directory does not exist: " << blast_dir << endl;
    char command[500];
    sprintf(command, "mkdir %s", blast_dir);
    if (system(command) != 0) {
      cout << "cannot create the specified blast directory " << blast_dir
           << endl;
      exit(0);
    }
  } else if (S_ISREG(buf.st_mode)) {
    cout << "Error: you have given a file name as blast directory" << endl;
    cout << "       please give a writable directory name " << blast_dir
         << endl;
    exit(0);
  }
  // check if the blast_dir exists, if it does not, create it
  else if (!S_ISDIR(buf.st_mode)) {
    char command[500];
    sprintf(command, "mkdir %s", blast_dir);
    if (system(command) != 0) {
      cout << "cannot create the specified blast directory " << blast_dir
           << endl;
      exit(0);
    }
  }
  // check if the directory has read and write permission
  if (access(blast_dir, R_OK) < 0) {
    cout << "Error: blast directory not readable " << blast_dir << endl;
    exit(0);
  }
  if (access(blast_dir, W_OK) < 0) {
    cout << "Error: blast directory not writable " << blast_dir << endl;
    exit(0);
  }
  if (access(blast_dir, X_OK) < 0) {
    cout << "Error: blast directory not executable " << blast_dir << endl;
    exit(0);
  }

  // cout << "Setting up a blast directory: " << blast_dir << endl << endl;
}

void printParameters() {
  cout << "List of parameters: " << endl;
  // cout << "  ave_grp_thr: " << ave_grp_thr << endl;
  // cout << "  minprob: " << minProb << endl;
  // cout << "  outfile name: " << outFile << endl;
  // cout << "  using probscons blosum62: " << probconsBLOSUM << endl;
  // cout << "  using local alignment option: " << useLocal;
  /*switch (useLocal) {
          case 1: cout << "  local only" << endl; break;
          case 0: cout << "  global only" << endl; break;
          case 2: cout << "  both local and global" << endl;
                  cout << "  weight of global: " << weightG << endl; break;
          default: cout << "   global only" << endl;
  }*/
  // cout << "  unaligned: " << unaligned << endl;
  // cout << "  ss: " << ss << "\n  solv: " << solv << "\n  parameter file: " <<
  // parameter_file << endl;
  cout << "  identity threshold: " << id_thr << endl;
  // cout << "  relax_number: " << relax_number << endl;
  // cout << "  reverse align order: " << reverse_align_order << endl;
}

void printParameters1() {
  // cout << "  psipred_dir: " << psipred_dir << endl;
  // cout << "  psipred_env_number:  " << psipred_env_number << endl;
  // cout << "  psipred_param: " << psipred_parameter_file << endl;
  cout << "  blast_dir: " << blast_dir << endl;
  // cout << "  use_ss_freq: " << use_ss_freq << endl;
  cout << "  secondary structure weight: " << ss_w << endl;
  cout << "  amino acid weight: " << score_w << endl;
  // cout << "  adjust_weight: " << adjust_weight << endl;
  // cout << "  score_shift: " << score_shift << endl;

  // cout << endl;
  // cout << "  blastpgp_cmd: " << blastpgp_command << endl;
  // cout << "  runpsipred_command: " << runpsipred_command << endl;
  // cout << "  runpsipred1_command: " << runpsipred1_command << endl;
  // cout << "  uniref90_file: " << uniref90_file << endl;
  // cout << endl;

  // cout << "  relax_count: ";
  // if(relax_count<=0) cout << "20" << endl;
  // else cout << relax_count << endl;
}

void printHelp(int prog) {
  if (prog == 1)  // mummals
  {
    cout << endl << " promals with 3D information " << endl;

    cout << " Usage:" << endl;
    cout << " \t promals input_file [options]" << endl;
    cout << endl;
    return;

    cout << " options: " << endl << endl;

    // cout << " -ss         Number of secondary structural types, 1 or 3" <<
    // endl; cout << " -solv       Number of solvent accessibility categories, 1,
    // 2 or 3" << endl; cout << " -unaligned  Number of additional match states
    // for unaligned regions, 0 or 1" << endl; cout << " -param      Input
    // parameter file for hidden Markov model" << endl;
    cout << " -id_thr     identity threshold above which neighboring groups "
            "are aligned in a fast way, between 0 and 1"
         << endl;

    cout << " -outfile    Output file name" << endl;
    cout << endl;

    cout << " Example: " << endl << endl;
    cout << " promals tmp1.fa -id_thr 0.6 -outfile tmp1.mummals.aln" << endl;

    cout << endl;
  }

  else if (prog == 0)  // merge_align
  {
    cout << endl
         << "merge_align - merge several multiple alignments into one "
            "alignment based on consistency";
    cout << endl << endl;

    cout << " input should be a file that contains a list of alignment file "
            "names"
         << endl;
    cout << " each alignment file should be in fasta format " << endl;
    cout << endl;
    cout << " Usage:" << endl;
    cout << " \t merge_align input_file_list [options]" << endl;
    cout << endl;

    cout << " options: " << endl << endl;
    cout << " -outfile    Output file name " << endl;
    cout << "	   if not specified, it will be input_file_list.merge_aln.aln"
         << endl;

    cout << endl;
  }
}
