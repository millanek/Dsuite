//
//  Dsuite_fBranch.cpp
//  DsuiteXcode
//
//  Created by Milan Malinsky on 11/11/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "Dsuite_fBranch.h"
#define SUBPROGRAM "Dbranch"

#define DEBUG 0

static const char *BRANCHSCORE_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] TREE_FILE.nwk DVALS_tree.txt\n"
"Implements the 'f-branch' type calculations developed by Hannes Svardal for Malinsky et al., 2018, Nat. Ecol. Evo.\n"
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
    std::cout << "Here we are fine" << std::endl;
    testTree->updateProgenyIds();
    std::cout << "Here we are still fine" << std::endl;
    for (std::vector<Branch*>::iterator b = testTree->branch.begin(); b != testTree->branch.end(); b++) {
        std::cout << (*b)->id << std::endl;
        std::cout << "(*b)->terminalSpeciesId: " << (*b)->terminalSpeciesId << std::endl;
        std::cout << "(*b)->daughterId1: " << (*b)->daughterId1 << std::endl;
        print_vector((*b)->progenyIds, std::cout);
        std::cout << std::endl;
    }
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


void Tree::updateProgenyIds() {
    // Determine the progeny of each branch (needed to know whether conditions are met, and for fossil constraints).
    // First of all, set progeniesComplete to 2 for all extinct and present branches.
    for (std::vector<Branch*>::iterator b = branch.begin(); b != branch.end(); b++) {
        if ((*b)->daughterId1 == "none") {
            (*b)->progeniesComplete = 2;
            (*b)->progenyIds.push_back((*b)->terminalSpeciesId);
        }
        // Set progenyPassedOn to true for the two root branches.
        if ((*b)->parentId == "treeOrigin") (*b)->progenyPassedOn = true;
    }
    bool allProgeniesComplete = false;
    while(!allProgeniesComplete) {
        std::vector<Branch*> newlyCompleted;
        for (std::vector<Branch*>::iterator b = branch.begin(); b != branch.end(); b++) {
            // Determine if the progeny of this branch is clear but has not been passed on to the parent yet.
            if ((*b)->progeniesComplete == 2 && (*b)->progenyPassedOn == false) {
                newlyCompleted.push_back(*b);
            }
        }
        if (newlyCompleted.size() == 0) allProgeniesComplete = true;
        for (std::vector<Branch*>::iterator b = newlyCompleted.begin(); b != newlyCompleted.end(); b++) {
            // Find parent, pass progeny+self on to parents progeny, add parent.progeniesComplete += 1, and change own progenyPassedOn to true.
            for (std::vector<Branch*>::iterator bb = branch.begin(); bb != branch.end(); bb++) {
                if ((*bb)->id == (*b)->parentId) {
                    (*bb)->progenyIds.insert((*bb)->progenyIds.end(), (*b)->progenyIds.begin(), (*b)->progenyIds.end() );
                    (*bb)->progeniesComplete++;
                    (*b)->progenyPassedOn = true;
                    break;
                }
            }
        }
    }
}
