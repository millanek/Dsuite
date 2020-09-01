//
//  Dsuite_common.h
//  DsuiteXcode
//
//  Created by Milan Malinsky on 21/07/2020.
//  Copyright © 2020 Milan Malinsky. All rights reserved.
//

#ifndef Dsuite_common_h
#define Dsuite_common_h

#include "Dsuite_utils.h"

void process_SETS_file(std::ifstream* setsFile, const string fName, std::map<string, std::vector<string>>& speciesToIDsMap, std::map<string, string>& IDsToSpeciesMap, int outgroupRequirement);
string makeHeader(bool quartet, bool includeFstats);
string prepareOutFileRootString(const string& providedPrefix, const string& runName, const string& setsFileName, const int regionStart, const int regionLength);
void printMissingLikelihoodsWarning(const string& chr, const string& pos);
void duplicateTreeValueError(const string& duplicate);
void assignTreeLevelsAndLinkToTaxa(string& treeLine, std::map<string,std::vector<int>>& taxaToLoc, std::vector<int>& levels);

#endif /* Dsuite_common_h */
