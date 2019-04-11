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

#define AUTHOR "Milan Malinsky"
#define PACKAGE_VERSION "0.1 r1"


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
"           abba-baba           Various forms of the D and f4 statistics for specific set of populations, also in windows along the genome\n"
"           Dmin                Calculate D statistics (ABBA-BABA) for all possible trios of samples\n"
"           DminCombine         Combine results from different Dmin runs (e.g. per-chromosome)\n"
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
        
        if(command == "abba-baba")
            abbaBabaMain(argc - 1, argv + 1);
        else if (command == "Dmin")
            DminMain(argc - 1, argv + 1);
        else if (command == "DminCombine")
            DminCombineMain(argc - 1, argv + 1);
        else
        {
            std::cerr << "Unrecognized command: " << command << "\n";
            return 1;
        }
        return 0;
    }
}

