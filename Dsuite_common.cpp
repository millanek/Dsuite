//
//  Dsuite_common.cpp
//  DsuiteXcode
//
//  Created by Milan Malinsky on 21/07/2020.
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

string makeHeader(bool quartet, bool includeFstats) {
    string header = "P1\tP2\tP3"; if (quartet) header += "\tP4";
    header += "\tDstatistic\tZ-score\tp-value"; if (includeFstats) { header += "\t"; header += F4HEADER; }
    header += "\tBBAA\tABBA\tBABA";
    return header;
}

string prepareOutFileRootString(const string& providedPrefix, const string& runName, const string& setsFileName, const int regionStart, const int regionLength) {
    string fileNameRootString; string outRoot; if (providedPrefix == "") { outRoot = stripExtension(setsFileName);} else { outRoot = providedPrefix; }
    if (regionStart == -1) { if (runName != "") fileNameRootString = outRoot + "_" + runName; else fileNameRootString = outRoot; }
    else fileNameRootString = outRoot+"_"+runName+"_"+numToString(regionStart)+"_"+numToString(regionStart+regionLength);
    return fileNameRootString;
}

void printMissingLikelihoodsWarning(const string& chr, const string& pos) {
    std::cerr << "WARNING: Could not fing genotype likelihoods/probabilities (GP, PL, or GL fields) for variant at " << chr << " " << pos << std::endl;
    std::cerr << "WARNING: Did you really mean to use the -g option? Reverting to using called genotypes." << std::endl;
}

void duplicateTreeValueError(const string& duplicate) {
    std::cerr << "ERROR: Duplicate value in the tree \"" << duplicate << "\"\n";
    std::cerr << "Exiting\n";
    exit(1);
}

void printInitialMessageTriosQuartets(const int regionLengthOpt, const int VCFlineCount, const int JKblockSizeBasedOnNum, const int jkWindowSizeOpt, const int jkNumOpt) {
    if (regionLengthOpt > 0) { std::cerr << "The VCF region to be analysed contains " << VCFlineCount << " variants\n"; }
    else { std::cerr << "The VCF contains " << VCFlineCount << " variants\n"; }
    if (jkWindowSizeOpt == 0) std::cerr << "Going to use block size of " << JKblockSizeBasedOnNum << " variants to get " << jkNumOpt << " Jackknife blocks\n";
}

void assignTreeLevelsAndLinkToTaxa(string& treeLine, std::map<string,std::vector<int>>& taxaToLoc, std::vector<int>& levels) {
    // First take care of any branch lengths
    std::regex branchLengths(":.*?(?=,|\\))");
    treeLine = std::regex_replace(treeLine,branchLengths,"");
    //std::cerr << line << std::endl;

    // Now process the tree
    levels.assign(treeLine.length(),0); int currentLevel = 0;
    std::vector<string> treeTaxonNames;
    string currentTaxonName = "";
    int lastBegin = 0;
    for (int i = 0; i < treeLine.length(); ++i) {
        if (treeLine[i] == '(') {
            currentLevel++; levels[i] = currentLevel;
        } else if (treeLine[i] == ')') {
            currentLevel--; levels[i] = currentLevel;
            if (currentTaxonName != "") {
                if (taxaToLoc.count(currentTaxonName) == 1) { duplicateTreeValueError(currentTaxonName); }
                treeTaxonNames.push_back(currentTaxonName);
                taxaToLoc[currentTaxonName].push_back(lastBegin);
                taxaToLoc[currentTaxonName].push_back(i-1);
                currentTaxonName = "";
            }
        } else if (treeLine[i] == ',') {
            levels[i] = currentLevel;
            if (currentTaxonName != "") {
                treeTaxonNames.push_back(currentTaxonName);
                taxaToLoc[currentTaxonName].push_back(lastBegin);
                taxaToLoc[currentTaxonName].push_back(i-1);
                currentTaxonName = "";
            }
        } else {
            if (currentTaxonName == "")
                lastBegin = i;
            levels[i] = currentLevel;
            currentTaxonName += treeLine[i];
        }
    }
    //print_vector(treeTaxonNames, std::cout,'\n');
    //print_vector(treeLevels, std::cout,' ');
    //for (std::map<string,std::vector<int>>::iterator i = treeTaxonNamesToLoc.begin(); i != treeTaxonNamesToLoc.end(); i++) {
    //    std::cout << i->first << "\t" << i->second[0] << "\t" << i->second[1] << "\t" << treeLevels[i->second[0]] << "\t" << treeLevels[i->second[1]] << std::endl;
    //}
}

int assignNumLinesToAnalyse(const int providedNumLinesOpt, const int regionLengthOpt,const string& vcfFileOpt) {
    int VCFlineCount;
    if (providedNumLinesOpt > 0) {
        VCFlineCount = providedNumLinesOpt;
    } else if (regionLengthOpt > 0) {
        VCFlineCount = regionLengthOpt;
    } else { // Block to find the number of lines in the VCF file
        std::istream* vcfFile = createReader(vcfFileOpt.c_str());
        // See how big is the VCF file
        vcfFile->unsetf(std::ios_base::skipws); // new lines will be skipped unless we stop it from happening:
        // count the newlines with an algorithm specialized for counting:
        VCFlineCount = (int)std::count(std::istream_iterator<char>(*vcfFile),std::istream_iterator<char>(),'\n');
        //std::cout << "VCF Lines: " << VCFlineCount << "\n";
    }
    return VCFlineCount;
}
