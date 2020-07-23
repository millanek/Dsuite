//
//  Dsuite_common.cpp
//  DsuiteXcode
//
//  Created by Milan Malinsky on 21/07/2020.
//  Copyright © 2020 Milan Malinsky. All rights reserved.
//

#include "Dsuite_common.h"

void process_SETS_file(std::ifstream* setsFile, const string fName, std::map<string, std::vector<string>>& speciesToIDsMap, std::map<string, string>& IDsToSpeciesMap, int outgroupRequirement) {
    int l = 0; string line;
    bool outgroupSpecified = false;
    while (getline(*setsFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        // std::cerr << line << std::endl;
        l++; if (line == "") { std::cerr << "Please fix the format of the " << fName << " file.\nLine " << l << " is empty." << std::endl; exit(EXIT_FAILURE); }
        std::vector<string> ID_Species = split(line, '\t');
        if (ID_Species.size() != 2) { std::cerr << "Please fix the format of the " << fName << " file.\nLine " << l << " does not have two columns separated by a tab." << std::endl; exit(EXIT_FAILURE); }
        if (ID_Species[1] == "Outgroup") { outgroupSpecified = true; }
        speciesToIDsMap[ID_Species[1]].push_back(ID_Species[0]);
        IDsToSpeciesMap[ID_Species[0]] = ID_Species[1];
        //std::cerr << ID_Species[1] << "\t" << ID_Species[0] << std::endl;
    }
    
    // Provide error/warning messages depending on which analysis is run and the presence/absence of Outgroup in the SETS file
    if (outgroupRequirement == OutgroupNotRequired && outgroupSpecified) { std::cerr << "WARNING: You specified the \"Outgroup\" in " << fName << ". This is needed in Dtrios, but will be ignored in Dquarters - the \"Outgroup\" will be treated as any other population. It must also be present in the tree if you are supplying one." << std::endl;
    }
    
    if (outgroupRequirement == OutgroupRequired && !outgroupSpecified) {
        std::cerr << "ERROR: The file " << fName << " needs to specify the \"Outgroup\"" << std::endl; exit(1);
    }
}

string makeHeader(bool includeFstats) {
    string header = "P1\tP2\tP3\tP4\tDstatistic\tZ-score\tp-value"; if (includeFstats) header += "\tf4-ratio";
    header += "\tBBAA\tABBA\tBABA";
    return header;
}

string prepareOutFileRootString(const string& providedPrefix, const string& runName, const string& setsFileName, const int regionStart, const int regionLength) {
    string fileNameRootString; string outRoot; if (providedPrefix == "") { outRoot = stripExtension(setsFileName);} else { outRoot = providedPrefix; }
    if (regionStart == -1) { if (runName != "") fileNameRootString = outRoot + "_" + runName; else fileNameRootString = outRoot; }
    else fileNameRootString = outRoot+"_"+runName+"_"+numToString(regionStart)+"_"+numToString(regionStart+regionLength);
    return fileNameRootString;
}
