//
//  Dsuite.cpp
//  Dsuite
//
//  Created by Milan Malinsky on 02/04/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include <iostream>
#include "Dsuite_utils.h"
#include "Dmin.h"
#include "D.h"
#include "Dmin_combine.h"
#include "Dsuite_fBranch.h"

#define AUTHOR "Milan Malinsky"
#define PACKAGE_VERSION "0.2 r19"


static const char *VERSION_MESSAGE =
"Dsuite software Version " PACKAGE_VERSION "\n"
"Written by Milan Malinsky.\n"
"\n";

static const char *USAGE_MESSAGE =
"Program: " PROGRAM_BIN "\n"
"Version: " PACKAGE_VERSION "\n"
"Contact: " AUTHOR " [" PACKAGE_BUGREPORT "]\n"
"Usage: " PROGRAM_BIN " <command> [options]\n\n"
"Commands:\n"
"           Dtrios                  Calculate D statistics (ABBA-BABA) for all possible trios of populations/species\n"
"           DtriosCombine           Combine results from Dtrios runs across genomic regions (e.g. per-chromosome)\n"
"           Dinvestigate            Follow up analyses for trios with significantly elevated D:\n"
"                                   calculates the f4 statistic, and also f_d and f_dM in windows along the genome\n"
"           Fbranch                 Calculate D and f statistics for branches on a tree that relates the populations/species\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

int main(int argc, char **argv) {
    
    if(argc <= 1)
    {
        std::cout << USAGE_MESSAGE;
        return 0;
    }
    else
    {
        std::string command(argv[1]);
        if(command == "help" || command == "--help" || command == "-h")
        {
            std::cout << USAGE_MESSAGE;
            return 0;
        }
        else if(command == "version" || command == "--version")
        {
            std::cout << VERSION_MESSAGE;
            return 0;
        }
        
        if(command == "Dinvestigate")
            abbaBabaMain(argc - 1, argv + 1);
        else if (command == "Dtrios")
            DminMain(argc - 1, argv + 1);
        else if (command == "DtriosCombine")
            DminCombineMain(argc - 1, argv + 1);
        else if (command == "Fbranch")
            fBranchMain(argc - 1, argv + 1);
        else
        {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
        return 0;
    }
}

