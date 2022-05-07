//
//  Dsuite_common.h
//  DsuiteXcode
//
//  Created by Milan Malinsky on 21/07/2020.
//

#ifndef Dsuite_common_h
#define Dsuite_common_h

#include "Dsuite_utils.h"

void process_SETS_file(std::ifstream* setsFile, const string fName, std::map<string, std::vector<string>>& speciesToIDsMap, std::map<string, string>& IDsToSpeciesMap, int outgroupRequirement);
string makeHeader(bool quartet, bool includeFstats, bool includeKSstats);
string prepareOutFileRootString(const string& providedPrefix, const string& runName, const string& setsFileName, const int regionStart, const int regionLength);
void printMissingLikelihoodsWarning(const string& chr, const string& pos);
void printInitialMessageTriosQuartets(const int regionLengthOpt, const int VCFlineCount, const int JKblockSizeBasedOnNum, const int jkWindowSizeOpt, const int jkNumOpt);
void duplicateTreeValueError(const string& duplicate);
void assignTreeLevelsAndLinkToTaxa(string& treeLine, std::map<string,std::vector<int>>& taxaToLoc, std::vector<int>& levels);
int assignNumLinesToAnalyse(const int providedNumLinesOpt, const int regionLengthOpt,const string& vcfFileOpt);

#endif /* Dsuite_common_h */
