//
//  Dmin_combine.cpp
//  Dsuite
//
//  Created by Milan Malinsky on 11/04/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "Dmin_combine.h"
#include "Dsuite_common.h"

#define SUBPROGRAM "DtriosCombine"

#define DEBUG 1

static const char *DMINCOMBINE_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] DminFile1 DminFile2 DminFile3 ....\n"
"Combine the BBAA, ABBA, and BABA counts from multiple files (e.g per-chromosome) and output the overall D stats,\n"
"p-values and f4-ratio values\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"       -t , --tree=TREE_FILE.nwk               (optional) a file with a tree in the newick format specifying the relationships between populations/species\n"
"                                               D and f4-ratio values for trios arranged according to the tree will be output in a file with _tree.txt suffix\n"
"       -s , --subset=start,length              (optional) only process a subset of the trios\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


static const char* shortopts = "hn:t:s:";

static const struct option longopts[] = {
    { "subset",   required_argument, NULL, 's' },
    { "run-name",   required_argument, NULL, 'n' },
    { "tree",   required_argument, NULL, 't' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static std::vector<string> dminFiles;
    static string runName = "";
    static string treeFile = "";
    int subsetStart = -1;
    int subsetLength = -1;
}


int DminCombineMain(int argc, char** argv) {
    parseDminCombineOptions(argc, argv);
    string line; // for reading the input files
    
    std::vector<std::istream*> dminstdErrFiles; std::vector<std::istream*> dminBBAAscoreFiles;
    for (int i = 0; i < opt::dminFiles.size(); i++) {
        std::istream* dminBBAAscoreFile;
        if (file_exists(opt::dminFiles[i] + "_combine.txt")) {
            dminBBAAscoreFile = createReader((opt::dminFiles[i] + "_combine.txt").c_str());
        } else if(file_exists(opt::dminFiles[i] + "_combine.txt.gz")) {
            dminBBAAscoreFile = createReader((opt::dminFiles[i] + "_combine.txt.gz").c_str());
        } else {
            std::cerr << "Can't find the file: " << opt::dminFiles[i] + "_combine.txt" << " or " << opt::dminFiles[i] + "_combine.txt.gz. Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
        dminBBAAscoreFiles.push_back(dminBBAAscoreFile);
        std::istream* dminstdErrFile;
        if (file_exists(opt::dminFiles[i] + "_combine_stderr.txt")) {
            dminstdErrFile = createReader((opt::dminFiles[i] + "_combine_stderr.txt").c_str());
        } else if(file_exists(opt::dminFiles[i] + "_combine_stderr.txt.gz")) {
            dminstdErrFile = createReader((opt::dminFiles[i] + "_combine_stderr.txt.gz").c_str());
        } else {
            std::cerr << "Can't find the file: " << opt::dminFiles[i] + "_combine_stderr.txt" << " or " << opt::dminFiles[i] + "_combine_stderr.txt.gz. Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
        dminstdErrFiles.push_back(dminstdErrFile);
        std::cerr << "Reading file " << opt::dminFiles[i] << std::endl;
    }
    
    std::istream* treeFile; std::ofstream* outFileTree;
    std::map<string,std::vector<int>> treeTaxonNamesToLoc; std::vector<int> treeLevels;
    if (opt::treeFile != "") {
        treeFile = new std::ifstream(opt::treeFile.c_str());
        if (!treeFile->good()) { std::cerr << "The file " << opt::treeFile << " could not be opened. Exiting..." << std::endl; exit(1);}
        outFileTree = new std::ofstream(opt::runName + "combined_tree.txt");
        getline(*treeFile, line);
        assignTreeLevelsAndLinkToTaxa(line,treeTaxonNamesToLoc,treeLevels);
    }
    // Now get the standard error values
    std::ofstream* outFileBBAA = new std::ofstream(opt::runName + "combined_BBAA.txt"); std::ofstream* outFileDmin = new std::ofstream(opt::runName + "combined_Dmin.txt");
    
    std::vector<double> BBAA_local_Ds; std::vector<double> ABBA_local_Ds; std::vector<double> BABA_local_Ds;
    string s1; string s2; string s3;
    bool allDone = false; bool fIncluded = false;
    int processedTriosNumber = 0; int exceptionCount = 0;
    
    getline(*dminBBAAscoreFiles[0], line); std::vector<string> patternCounts = split(line, '\t');
    if (patternCounts.size() == 12) fIncluded = true;
    string header = makeHeader(fIncluded);
    *outFileBBAA << header << std::endl; *outFileDmin << header << std::endl;
    if (opt::treeFile != "") *outFileTree << header << std::endl;
    dminBBAAscoreFiles[0]->seekg(0, dminBBAAscoreFiles[0]->beg); // Go back to the beginning of this file
    
    do {
        TrioDinfo info; processedTriosNumber++;
        if (processedTriosNumber % 10000 == 0) { std::cerr << "Processed " << processedTriosNumber << " trios" << std::endl; }
        
        if (opt::subsetStart != -1) {
            if (processedTriosNumber < opt::subsetStart) {
                for (int i = 0; i < dminBBAAscoreFiles.size(); i++) { getline(*dminBBAAscoreFiles[i], line); }
                for (int i = 0; i < dminstdErrFiles.size(); i++) { getline(*dminstdErrFiles[i], line); }
                continue;
            }
            if (processedTriosNumber >= (opt::subsetStart+opt::subsetLength)) {
                std::cerr << "DONE" << std::endl; break;
            }
        }
        
        
        for (int i = 0; i < dminBBAAscoreFiles.size(); i++) {
            if (getline(*dminBBAAscoreFiles[i], line)) {
                std::vector<string> patternCounts = split(line, '\t');
                assert(patternCounts.size() == 6 || patternCounts.size() == 12);

                if (i == 0) {
                    s1 = patternCounts[0]; s2 = patternCounts[1]; s3 = patternCounts[2];
                } else {
                    assert(s1 == patternCounts[0]); assert(s2 == patternCounts[1]); assert(s3 == patternCounts[2]);
                }
                info.BBAAtotal += stringToDouble(patternCounts[3]);
                info.BABAtotal += stringToDouble(patternCounts[4]);
                info.ABBAtotal += stringToDouble(patternCounts[5]);
                if (fIncluded) {
                    info.F_G_denom1 += stringToDouble(patternCounts[6]);
                    info.F_G_denom2 += stringToDouble(patternCounts[7]);
                    info.F_G_denom3 += stringToDouble(patternCounts[8]);
                    info.F_G_denom1_reversed += stringToDouble(patternCounts[9]);
                    info.F_G_denom2_reversed += stringToDouble(patternCounts[10]);
                    info.F_G_denom3_reversed += stringToDouble(patternCounts[11]);
                }
            } else {
                allDone = true; break;
            }
        }
        
        for (int i = 0; i < dminstdErrFiles.size(); i++) {
            if (getline(*dminstdErrFiles[i], line)) {
                std::vector<string> localDs = split2(line, "\t");
                //assert(localDs.size() == 3 || localDs.size() == 0);
                if (localDs.size() == 3) {
                    std::vector<string> regionD_strings0 = split(localDs[0], ',');
                    std::vector<string> regionD_strings1 = split(localDs[1], ',');
                    std::vector<string> regionD_strings2 = split(localDs[2], ',');
                    for (int j = 0; j < regionD_strings0.size(); j++) {
                        double localD = stringToDouble(regionD_strings0[j]);
                        if (!std::isnan(localD)) info.regionDs[0].push_back(localD);
                    }
                    for (int j = 0; j < regionD_strings1.size(); j++) {
                        double localD = stringToDouble(regionD_strings1[j]);
                        if (!std::isnan(localD)) info.regionDs[1].push_back(localD);
                    }
                    for (int j = 0; j < regionD_strings2.size(); j++) {
                        double localD = stringToDouble(regionD_strings2[j]);
                        if (!std::isnan(localD)) info.regionDs[2].push_back(localD);
                    }
                } else {
                    print_vector(localDs,std::cerr); exit(EXIT_FAILURE);
                }
            } else {
                allDone = true; break;
            }
        }
        
        
        if (!allDone) {
            try {
                info.calculateFinalDs();
            } catch (const char* msg) {
                exceptionCount++;
                if (exceptionCount <= 10) {
                    std::cerr << msg << std::endl;
                    std::cerr << "Could not calculate p-values for the trio: " << s1 << " " << s2 << " " << s3 << std::endl;
                    std::cerr << "You should probably decrease the the jackknife block size (-j option)" << std::endl;
                    std::cerr << std::endl;
                }
                info.D1_p = nan(""); info.D2_p = nan(""); info.D3_p = nan("");
            }
            
            std::vector<string> trio; trio.push_back(s1); trio.push_back(s2); trio.push_back(s3);
            // Find which topology is in agreement with the counts of BBAA, BABA, and ABBA
            info.assignBBAAarrangement();
            std::vector<string> BBAAoutVec = info.makeOutVec(trio, fIncluded, info.BBAAarrangement);
            print_vector(BBAAoutVec,*outFileBBAA);
           
            // Find Dmin:
            info.assignDminArrangement();
            std::vector<string> DminOutVec = info.makeOutVec(trio, fIncluded, info.DminArrangement);
            print_vector(DminOutVec,*outFileDmin);
            
            if (opt::treeFile != "") {
                int loc1 = treeTaxonNamesToLoc[s1][0]; int loc2 = treeTaxonNamesToLoc[s2][0]; int loc3 = treeTaxonNamesToLoc[s3][0];
                info.treeArrangement = info.assignTreeArrangement(treeLevels, loc1, loc2, loc3);
                std::vector<string> treeOutVec = info.makeOutVec(trio, fIncluded, info.treeArrangement);
                print_vector(treeOutVec,*outFileTree);
            }
        }
 
    } while(!allDone);
    
    if (exceptionCount > 10) {
        std::cerr << "..." << std::endl;
        std::cerr << "p-value could not be calculated for " << exceptionCount << " trios" << std::endl;
        std::cerr << "You should definitely decrease the the jackknife block size!!!" << std::endl;
        std::cerr << std::endl;
    }
    
    return 0;
    
}



void parseDminCombineOptions(int argc, char** argv) {
    bool die = false; string subsetArgString; std::vector<string> subsetArgs;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 't': arg >> opt::treeFile; break;
            case 's': arg >> subsetArgString; subsetArgs = split(subsetArgString, ',');
                opt::subsetStart = (int)stringToDouble(subsetArgs[0]); opt::subsetLength = (int)stringToDouble(subsetArgs[1]);  break;
            case 'h':
                std::cout << DMINCOMBINE_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    
    int nFilenames = argc - optind;
    if (nFilenames < 1) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DMINCOMBINE_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    while (optind < argc) {
        opt::dminFiles.push_back(argv[optind++]);
    }
}
