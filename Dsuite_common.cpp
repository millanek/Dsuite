//
//  Dsuite_common.cpp
//  DsuiteXcode
//
//  Created by Milan Malinsky on 21/07/2020.
//

#include "Dsuite_common.h"



void SetInformation::linkSetsAndVCFpositions(const std::vector<std::string>& sampleNames) {
    // print_vector_stream(sampleNames, std::cerr);
    for (std::vector<std::string>::size_type i = 0; i != sampleNames.size(); i++) {
        try { posToPopMap[i] = IDsToPopMap.at(sampleNames[i]); } catch (const std::out_of_range& oor) {
            std::cerr << "WARNING: The sample " << sampleNames[i] << " is in the VCF but not assigned in the SETS.txt file" << std::endl;
        }
    }
    // Iterate over all the keys in the map to find the samples in the VCF:
    // Give an error if no sample is found for a species:
    for(std::map<string, std::vector<string>>::const_iterator it = popToIDsMap.begin(); it != popToIDsMap.end(); ++it) {
        string sp =  it->first;
        //std::cerr << "sp " << sp << std::endl;
        std::vector<string> IDs = it->second;
        std::vector<size_t> spPos = locateSet(sampleNames, IDs); 
        if (spPos.empty()) {
            std::cerr << "Did not find any samples in the VCF for \"" << sp << "\"" << std::endl;
            assert(!spPos.empty());
        }
        popToPosMap[sp] = spPos;
    }
}


string makeHeader(bool quartet, bool includeFstats, bool includeKSstats) {
    string header = "P1\tP2\tP3"; if (quartet) header += "\tP4";
    header += "\tDstatistic\tZ-score\tp-value";
    if (includeFstats) { header += "\t"; header += F4HEADER; }
    if (includeKSstats) { header += "\t"; header += "clustering_KS_p-val1"; header += "\t"; header += "clustering_KS_p-val2";}
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
