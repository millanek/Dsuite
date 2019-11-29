//
//  D.h
//  Dsuite
//
//  Created by Milan Malinsky on 11/04/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#ifndef D_h
#define D_h

#include "Dsuite_utils.h"

class TestTrioInfo {
public:
    TestTrioInfo(int windowSize) {
        windowABBAs.resize(windowSize); windowBABAs.resize(windowSize);
        windowF_d_denoms.resize(windowSize); windowF_dM_denoms.resize(windowSize);
        windowInformativeSitesCords.resize(windowSize);
        interimF_d_denom = 0; interimF_dM_denom = 0;

        usedVars = 0;
        ABBAtotal = 0; BABAtotal = 0;
        F_d_denom = 0; F_dM_denom = 0;
        F_G_denom = 0; F_G_num = 0;

    };
    
    // string P1; string P2; string P3;
    
    std::deque<double> windowABBAs; std::deque<double> windowBABAs;
    std::deque<double> windowF_d_denoms; std::deque<double> windowF_dM_denoms;
    std::deque<int> windowInformativeSitesCords;
    double interimF_d_denom; double interimF_dM_denom;
    //double D1; double D2; double D3; double D1_p; double D2_p; double D3_p;
    
    double ABBAtotal; double BABAtotal;
    double F_d_denom; double F_dM_denom; double F_G_denom; double F_G_num;
    int usedVars;
    
    
    
};


void parseAbbaBabaOptions(int argc, char** argv);
int abbaBabaMain(int argc, char** argv);
#endif /* D_h */
