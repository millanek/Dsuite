//
//  Dsuite_common.h
//  DsuiteXcode
//
//  Created by Milan Malinsky on 21/07/2020.
//

#ifndef Dsuite_common_h
#define Dsuite_common_h

#define stdInInfo   "Use 'stdin' for the VCF file when piping from another program into Dsuite via standard input\n" \
                    "in this case it is necessary to provide the number of lines in the filtered VCF via the -l option\n" \
                    "For example, to filter the VCF for overall mimimum depth of at least 1000 across all samples:\n" \
                    "NUMLINES=$(bcftools view -i 'INFO/DP>1000' INPUT_FILE.vcf | wc -l)  # to get NUMLINES\n" \
                    "bcftools view -i 'INFO/DP>1000' INPUT_FILE.vcf | Dsuite Dtrios -l $NUMLINES stdin SETS.txt\n" \
                    "\n"

#define regionOption    "       -r, --region=start,length               (optional) only process a subset of the VCF file; both \"start\" and \"length\" indicate variant numbers\n" \
                        "                                               e.g. --region=20001,10000 will process variants from 20001 to 30000\n"

#define treeOption      "       -t, --tree=TREE_FILE.nwk                (optional) a file with a tree in the newick format specifying the relationships between populations/species\n" \
                        "                                               D and f4-ratio values for trios arranged according to the tree will be output in a file with _tree.txt suffix\n"

#define outOption       "       -o, --out-prefix=OUT_FILE_PREFIX        (optional) the prefix for the files where the results should be written\n" \
                        "                                               output will be put in OUT_FILE_PREFIX_BBAA.txt, OUT_FILE_PREFIX_Dmin.txt, OUT_FILE_PREFIX_tree.txt etc.\n" \
                        "                                               by default, the prefix is taken from the name of the SETS.txt file\n"

#include "Dsuite_utils.h"

inline void notEnoughPopulationsError(const int minPopulations) {
    std::cerr << "ERROR: You need at least " << minPopulations << " sets (populations/species) for this analysis." << std::endl;
    exit(EXIT_FAILURE);
}

inline void outgroupNeededError(const string& setsFileName) {
    std::cerr << "ERROR: The file " << setsFileName << " needs to specify the \"Outgroup\"" << std::endl;
    exit(EXIT_FAILURE);
}

inline void outgroupNotUsedInQuartetsWarning(const string& setsFileName) {
    std::cerr << "WARNING: You specified the \"Outgroup\" in " << setsFileName << ". This is needed in Dtrios, but will be ignored in Dquarters - the \"Outgroup\" will be treated as any other population. It must also be present in the tree if you are supplying one." << std::endl;
}

inline void wrongNumberOfColumnsError(const string& setsFileName, int lineNum) {
    std::cerr << "ERROR: Please fix the format of the " << setsFileName << " file." << std::endl;
    std::cerr << "Line " << lineNum << " does not have two columns separated by a tab." << std::endl;
    exit(EXIT_FAILURE);
}

inline void lineEmptyError(const string& setsFileName, int lineNum) {
    std::cerr << "ERROR: Please fix the format of the " << setsFileName << " file." << std::endl;
    std::cerr << "Line " << lineNum << " is empty." << std::endl;
    exit(EXIT_FAILURE);
}

class SetInformation {
public:
    
    SetInformation(const string& setsFileName, const int minPopulations, const int outgroupRequirement) {
        
        std::ifstream* setsFile = new std::ifstream(setsFileName.c_str());
        assertFileOpen(*setsFile, setsFileName);
        
        string line; int l = 0; bool outgroupSpecified = false;
        while (getline(*setsFile, line)) {
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
            
            l++; if (line == "") lineEmptyError(setsFileName,l);
            
            std::vector<string> ID_Pop = split(line, '\t');
            
            if (ID_Pop.size() != 2) wrongNumberOfColumnsError(setsFileName,l);
            if (ID_Pop[1] == "Outgroup") { outgroupSpecified = true; }
            
            popToIDsMap[ID_Pop[1]].push_back(ID_Pop[0]);
            IDsToPopMap[ID_Pop[0]] = ID_Pop[1];
        }
        
        for(std::map<string,std::vector<string>>::iterator it = popToIDsMap.begin(); it != popToIDsMap.end(); ++it) {
            if ((it->first) != "Outgroup" && it->first != "xxx") {
                populations.push_back(it->first);
            }
        } std::cout << "There are " << populations.size() << " sets (populations/species) excluding the Outgroup" << std::endl;
        
        if (populations.size() < minPopulations) notEnoughPopulationsError(minPopulations);
        
        // Provide error/warning messages depending on which analysis is run and the presence/absence of Outgroup in the SETS file
        if (outgroupRequirement == OutgroupNotRequired && outgroupSpecified) outgroupNotUsedInQuartetsWarning(setsFileName);
        if (outgroupRequirement == OutgroupRequired && !outgroupSpecified) outgroupNeededError(setsFileName);
    };
    
    
    
    std::vector<string> populations;
    std::map<string, string> IDsToPopMap;
    std::map<string, std::vector<string>> popToIDsMap;
    std::map<string, std::vector<size_t>> popToPosMap;
    std::map<size_t, string> posToPopMap;

    void linkSetsAndVCFpositions(const std::vector<std::string>& sampleNames);

};


void process_SETS_file(std::ifstream* setsFile, const string fName, std::map<string, std::vector<string>>& speciesToIDsMap, std::map<string, string>& IDsToSpeciesMap, int outgroupRequirement);
string makeHeader(bool quartet, bool includeFstats, bool includeKSstats);
string prepareOutFileRootString(const string& providedPrefix, const string& runName, const string& setsFileName, const int regionStart, const int regionLength);
void printMissingLikelihoodsWarning(const string& chr, const string& pos);
void printInitialMessageTriosQuartets(const int regionLengthOpt, const int VCFlineCount, const int JKblockSizeBasedOnNum, const int jkWindowSizeOpt, const int jkNumOpt);
void duplicateTreeValueError(const string& duplicate);
void assignTreeLevelsAndLinkToTaxa(string& treeLine, std::map<string,std::vector<int>>& taxaToLoc, std::vector<int>& levels);
int assignNumLinesToAnalyse(const int providedNumLinesOpt, const int regionLengthOpt,const string& vcfFileOpt);

inline void reportProgessVCF(const int variantsProcessed, const std::clock_t startTime) {
    double durationOverall = ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    std::cout << "Processed " << variantsProcessed << " variants in " << durationOverall << "secs" << std::endl;
}

inline void reportProgessVCF(const int variantsProcessed, const int VCFlineCount, const std::clock_t startTime) {
    double durationOverall = ( std::clock() - startTime ) / (double) CLOCKS_PER_SEC;
    std::cerr << "Processed " << variantsProcessed << " variants (" << ((double)variantsProcessed/VCFlineCount)*100 << "%) in " << durationOverall << "secs" << std::endl;
}

#endif /* Dsuite_common_h */
