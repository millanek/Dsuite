//
//  Dsuite_fBranch.cpp
//  DsuiteXcode
//
//  Created by Milan Malinsky on 11/11/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "Dsuite_fBranch.h"
#define SUBPROGRAM "Fbranch"

#define DEBUG 0

static const char *BRANCHSCORE_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] TREE_FILE.nwk FVALS_tree.txt\n"
"Implements the 'f-branch' type calculations developed by Hannes Svardal for Malinsky et al., 2018, Nat. Ecol. Evo.\n"
"Uses the f_G values produced by Dsuite Dtrios with the -f and --tree options; this is the output of Dtrios with the \"_tree.txt\" suffix\n"
"\n"
"       -h, --help                              display this help and exit\n"
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
}


int fBranchMain(int argc, char** argv) {
    parseFbranchOptions(argc, argv);
    std::istream* treeFile = new std::ifstream(opt::treeFile.c_str());
    if (!treeFile->good()) { std::cerr << "The file " << opt::treeFile << " could not be opened. Exiting..." << std::endl; exit(EXIT_FAILURE);}
    std::istream* DvalsFile = new std::ifstream(opt::DvalsFile.c_str());
    if (!DvalsFile->good()) { std::cerr << "The file " << opt::DvalsFile << " could not be opened. Exiting..." << std::endl; exit(EXIT_FAILURE);}
    if (opt::DvalsFile.substr(opt::DvalsFile.size()-9) != "_tree.txt") { std::cerr << "The name of the input file with the f_G values should end in \"_tree.txt\".\nPlease make sure you run Dtrios with the -f and --tree options and then feed the correct file into Fbranch. Exiting..." << std::endl; exit(EXIT_FAILURE); }
    std::map<string,std::vector<std::vector<string>>> acToBmap;
    string line; int l = 0;
    getline(*DvalsFile, line); // skip over the header
    while (getline(*DvalsFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        l++; if (line == "") { std::cerr << "Please fix the format of the " << opt::DvalsFile << " file.\nLine " << l << " is empty. Exiting..." << std::endl; exit(EXIT_FAILURE); }
        std::vector<string> speciesAndVals = split(line, '\t');
        if (speciesAndVals.size() <= 6) { std::cerr << "Please fix the format of the " << opt::DvalsFile << " file." << std::endl;
            std::cerr << "Looks like the file does not contain f statistics. Run Dsuite Dtrios with the -f option. Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
        std::vector<string> bAndVal;  bAndVal.push_back(speciesAndVals[1]); bAndVal.push_back(speciesAndVals[5]);
        std::vector<string> aAndVal;  aAndVal.push_back(speciesAndVals[0]); aAndVal.push_back(speciesAndVals[5]);
        acToBmap[speciesAndVals[0]+","+speciesAndVals[2]].push_back(bAndVal);
    }
    string treeString; getline(*treeFile, treeString);
    Tree* testTree = new Tree(treeString);
    testTree->updateProgenyIds();
    testTree->fillSisterBranches();
    for (std::vector<Branch*>::iterator b = testTree->branches.begin(); b != testTree->branches.end(); b++) {
        if ((*b)->parentId != "treeOrigin") {
            std::vector<string> Bs = (*b)->progenyIds;
            std::vector<string> As = (*b)->sisterBranch->progenyIds;
            std::vector<double> Bmins; std::vector<double> vals;
            for (std::vector<string>::iterator C = testTree->allSpecies.begin(); C != testTree->allSpecies.end(); C++) {
                for (std::vector<string>::iterator A = As.begin(); A != As.end(); A++) {
                    std::vector<std::vector<string>> bAndVal; std::vector<std::vector<string>> aAndVal;
                    try { bAndVal = acToBmap.at(*A+","+*C); } catch (const std::out_of_range& oor) {}
                    for (int i = 0; i < bAndVal.size(); i++) {
                        if (std::count(Bs.begin(), Bs.end(), bAndVal[i][0])) {
                            vals.push_back(stringToDouble(bAndVal[i][1]));
                        }
                    }
                    if (!vals.empty()) { Bmins.push_back(*std::min_element(vals.begin(),vals.end())); vals.clear(); }
 
                }
                double fbC = NAN;
                if (!Bmins.empty()) { fbC = median(Bmins.begin(),Bmins.end()); Bmins.clear(); }
                else { // There is no positive value; just find if any value is possible for this ABC combination
                    bool ACpossible = false;
                    for (std::vector<string>::iterator B = Bs.begin(); B != Bs.end(); B++) {
                        std::vector<std::vector<string>> bAndVal; std::vector<std::vector<string>> aAndVal;
                        try { bAndVal = acToBmap.at(*B+","+*C); } catch (const std::out_of_range& oor) {}
                        for (int i = 0; i < bAndVal.size(); i++) {
                            if (std::count(As.begin(), As.end(), bAndVal[i][0])) {
                                ACpossible = true; break;
                            }
                        }
                    }
                    if (ACpossible) fbC = 0;
                }
                (*b)->fbCvals.push_back(fbC);
            }
        }
    }
    std::cout << "branch\tbranch_descendants\t"; print_vector(testTree->allSpecies, std::cout);
    for (std::vector<Branch*>::iterator b = testTree->branches.begin(); b != testTree->branches.end(); b++) {
        if ((*b)->parentId != "treeOrigin") {
            std::cout << (*b)->id << "\t"; print_vector((*b)->progenyIds, std::cout, ',', false);
            std::cout << "\t"; print_vector((*b)->fbCvals, std::cout);
            //std::cout << "Sister branch:\t" <<  (*b)->sisterBranch->id << std::endl;
            //std::cout << "This branch progeny:\t"; print_vector((*b)->progenyIds, std::cout);
            //std::cout << "Sister branch progeny:\t"; print_vector((*b)->sisterBranch->progenyIds, std::cout);
            //std::cout << "fbCs:\t"; print_vector((*b)->fbCvals, std::cout);
            //std::cout << std::endl;
        }
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
    for (std::vector<Branch*>::iterator b = branches.begin(); b != branches.end(); b++) {
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
        for (std::vector<Branch*>::iterator b = branches.begin(); b != branches.end(); b++) {
            // Determine if the progeny of this branch is clear but has not been passed on to the parent yet.
            if ((*b)->progeniesComplete == 2 && (*b)->progenyPassedOn == false) {
                newlyCompleted.push_back(*b);
            }
        }
        if (newlyCompleted.size() == 0) allProgeniesComplete = true;
        for (std::vector<Branch*>::iterator b = newlyCompleted.begin(); b != newlyCompleted.end(); b++) {
            // Find parent, pass progeny+self on to parents progeny, add parent.progeniesComplete += 1, and change own progenyPassedOn to true.
            for (std::vector<Branch*>::iterator bb = branches.begin(); bb != branches.end(); bb++) {
                if ((*bb)->id == (*b)->parentId) {
                    (*b)->parentBranch = *bb;
                    (*bb)->progenyIds.insert((*bb)->progenyIds.end(), (*b)->progenyIds.begin(), (*b)->progenyIds.end() );
                    (*bb)->progeniesComplete++;
                    (*b)->progenyPassedOn = true;
                    break;
                }
            }
        }
    }
}

void Tree::fillSisterBranches() {
    for (std::vector<Branch*>::iterator b = branches.begin(); b != branches.end(); b++) {
        if ((*b)->parentId != "treeOrigin") {
            string sisterId;
            if ((*b)->parentBranch->daughterId1 != (*b)->id)
                sisterId = (*b)->parentBranch->daughterId1;
            else
                sisterId = (*b)->parentBranch->daughterId2;
            for (std::vector<Branch*>::iterator bb = branches.begin(); bb != branches.end(); bb++) {
                if ((*bb)->id == sisterId) {
                    (*b)->sisterBranch = *bb;
                    break;
                }
            }
        }
    }
}
