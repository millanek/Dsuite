//
//  Dsuite_fBranch.cpp
//  DsuiteXcode
//
//  Created by Milan Malinsky on 11/11/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "Dsuite_fBranch.h"
#define SUBPROGRAM "DbranchScore"

#define DEBUG 0

static const char *BRANCHSCORE_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] TREE_FILE.nwk DVALS_tree.txt\n"

"The test_trios.txt should contain names of three populations for which the statistics will be calculated:\n"
"POP1   POP2    POP3\n"
"There can be multiple lines and then the program generates multiple ouput files, named like POP1_POP2_POP3_localFstats_SIZE_STEP.txt\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -w SIZE,STEP --window=SIZE,STEP         (required) D, f_D, and f_dM statistics for windows containing SIZE useable SNPs, moving by STEP (default: 50,25)\n"
//"       --fJackKnife=WINDOW                     (optional) Calculate jackknife for the f_G statistic from Green et al. Also outputs \n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


//enum { OPT_F_JK };

static const char* shortopts = "hw:n:";

//static const int JK_WINDOW = 5000;

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "window",   required_argument, NULL, 'w' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string treeFile;
    static string DvalsFile;
    //int jkWindowSize = JK_WINDOW;
}


int fBranchMain(int argc, char** argv) {
    parseFbranchOptions(argc, argv);
    std::istream* treeFile = new std::ifstream(opt::treeFile.c_str());
    string treeString; getline(*treeFile, treeString);
    Tree* testTree = new Tree(treeString);
    return 0;
    
}

void parseFbranchOptions(int argc, char** argv) {
    bool die = false;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'h':
                std::cout << BRANCHSCORE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << BRANCHSCORE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::treeFile = argv[optind++];
    opt::DvalsFile = argv[optind++];
}
