//
//  Dmin.h
//  Dsuite
//
//  Created by Milan Malinsky on 02/04/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#ifndef Dmin_h
#define Dmin_h

#include "Dsuite_utils.h"
void parseDminOptions(int argc, char** argv);
int DminMain(int argc, char** argv);

class TrioDinfo {
public:
    TrioDinfo() {
        ABBAtotal = 0;
        BABAtotal = 0;
        BBAAtotal = 0;
        treeArrangement = 0; BBAAarrangement = 0; DminArrangement = 0;
        regionDs.resize(3);
    };
    
    // string P1; string P2; string P3;
    double ABBAtotal; double BABAtotal; double BBAAtotal;
    double D1; double D2; double D3; double D1_p; double D2_p; double D3_p;
    double F_d_denom; double F_d_denom_reversed;
    double F_dM_denom; double F_dM_denom_reversed;
    double F_G_denom; double F_G_denom_reversed;
    
    double localABBAtotal; double localBABAtotal; double localBBAAtotal; std::vector<std::vector<double>> regionDs; // vector with three empty (double) vectors
    
    int treeArrangement;    // 1 - trios[i][0] and trios[i][1] are P1 and P2
                            // 2 - trios[i][0] and trios[i][2] are P1 and P2
                            // 3 - trios[i][1] and trios[i][2] are P1 and P2
    
    int BBAAarrangement;    // 1 - trios[i][0] and trios[i][1] are P1 and P2
                            // 2 - trios[i][0] and trios[i][2] are P1 and P2
                            // 3 - trios[i][1] and trios[i][2] are P1 and P2
    
    int DminArrangement;    // 1 - trios[i][0] and trios[i][1] are P1 and P2
                            // 2 - trios[i][0] and trios[i][2] are P1 and P2
                            // 3 - trios[i][1] and trios[i][2] are P1 and P2
    
    
    void assignTreeArrangement(std::vector<int>& treeLevels, int loc1, int loc2, int loc3) {
        int midLoc = std::max(std::min(loc1,loc2), std::min(std::max(loc1,loc2),loc3));
        if (midLoc == loc1) {
            if (loc2 < loc1) {
                int m1 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc1);
                int m2 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc3);
                if (m1 < m2) treeArrangement = 2; else treeArrangement = 1;
            } else {
                int m1 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc1);
                int m2 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc2);
                if (m1 < m2) treeArrangement = 1; else treeArrangement = 2;
            }
        } else if (midLoc == loc2) {
            if (loc1 < loc2) {
                int m1 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc2);
                int m2 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc3);
                if (m1 < m2) treeArrangement = 3; else treeArrangement = 1;
            } else {
                int m1 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc2);
                int m2 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc1);
                if (m1 < m2) treeArrangement = 1; else treeArrangement = 3;
            }
        } else if (midLoc == loc3) {
            if (loc1 < loc3) {
                int m1 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc3);
                int m2 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc2);
                if (m1 < m2) treeArrangement = 3; else treeArrangement = 2;
                
            } else {
                int m1 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc3);
                int m2 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc1);
                if (m1 < m2) treeArrangement = 2; else treeArrangement = 3;
            }
        }
    }
    
    std::vector<string> assignBBAAarrangement(std::vector<string>& trio, bool fStats) {
        std::vector<string> outVec; if (fStats) outVec.resize(11); else outVec.resize(5);
        // Find which topology is in agreement with the counts of the BBAA, BABA, and ABBA patterns
        if (BBAAtotal >= BABAtotal && BBAAtotal >= ABBAtotal) {
            BBAAarrangement = 1; outVec[2] = trio[2];
            if (D1 >= 0) {outVec[0] = trio[0]; outVec[1] = trio[1]; }
            else { outVec[0] = trio[1]; outVec[1] = trio[0]; }
            outVec[3] = numToString(std::fabs(D1)); outVec[4] = numToString(D1_p);
        } else if (BABAtotal >= BBAAtotal && BABAtotal >= ABBAtotal) {
            BBAAarrangement = 2; outVec[2] = trio[1];
            if (D2 >= 0) {outVec[0] = trio[0]; outVec[1] = trio[2]; }
            else { outVec[0] = trio[2]; outVec[1] = trio[0]; }
            outVec[3] = numToString(std::fabs(D2)); outVec[4] = numToString(D2_p);
        } else if (ABBAtotal >= BBAAtotal && ABBAtotal >= BABAtotal) {
            BBAAarrangement = 3; outVec[2] = trio[0];
            if (D3 >= 0) {outVec[0] = trio[2]; outVec[1] = trio[1]; }
            else { outVec[0] = trio[1]; outVec[1] = trio[2]; }
            outVec[3] = numToString(std::fabs(D3)); outVec[4] = numToString(D3_p);
        }
        return outVec;
    }
    
    std::vector<string> assignDminArrangement(std::vector<string>& trio, bool fStats) {
        std::vector<string> outVec; if (fStats) outVec.resize(11); else outVec.resize(5);
        if (std::fabs(D1) <= std::fabs(D2) && std::fabs(D1) <= std::fabs(D3)) { // (P3 == S3)
            DminArrangement = 1; outVec[2] = trio[2];
            if (D1 >= 0) { outVec[0] = trio[0]; outVec[1] = trio[1]; }
            else {outVec[0] = trio[1]; outVec[1] = trio[0];}
            outVec[3] = numToString(std::fabs(D1)); outVec[4] = numToString(D1_p);
        } else if (std::fabs(D2) <= std::fabs(D1) && std::fabs(D2) <= std::fabs(D3)) { // (P3 == S2)
            DminArrangement = 2; outVec[2] = trio[1];
            if (D2 >= 0) {outVec[0] = trio[0]; outVec[1] = trio[2]; }
            else { outVec[0] = trio[2]; outVec[1] = trio[0]; }
            outVec[3] = numToString(std::fabs(D2)); outVec[4] = numToString(D2_p);
        } else if (std::fabs(D3) <= std::fabs(D1) && std::fabs(D3) <= std::fabs(D2)) { // (P3 == S1)
            DminArrangement = 3; outVec[2] = trio[0];
            if (D3 >= 0) {outVec[0] = trio[2]; outVec[1] = trio[1]; }
            else { outVec[0] = trio[1]; outVec[1] = trio[2]; }
            outVec[3] = numToString(std::fabs(D3)); outVec[4] = numToString(D3_p);
        }
        return outVec;
    }
    
    
    
    void calculateFinalDs() throw(const char*) {
        double* Ds = calculateThreeDs(ABBAtotal, BABAtotal, BBAAtotal);
        D1 = Ds[0]; D2 = Ds[1]; D3 = Ds[2];
        
        // Get the standard error values:
        double D1stdErr = jackknive_std_err(regionDs[0]); double D2stdErr = jackknive_std_err(regionDs[1]);
        double D3stdErr = jackknive_std_err(regionDs[2]);
        // Get the Z-scores
        double D1_Z = std::fabs(D1)/D1stdErr; double D2_Z = std::fabs(D2)/D2stdErr;
        double D3_Z = std::fabs(D3)/D3stdErr;
        // And p-values
        D1_p = 1 - normalCDF(D1_Z); D2_p = 1 - normalCDF(D2_Z);
        D3_p = 1 - normalCDF(D3_Z);
    }
    
    void addRegionDs() {
        double* Ds = calculateThreeDs(localABBAtotal, localBABAtotal, localBBAAtotal);
        regionDs[0].push_back(Ds[0]); regionDs[1].push_back(Ds[1]); regionDs[2].push_back(Ds[2]);
        localABBAtotal = 0; localBABAtotal = 0; localBBAAtotal = 0;
    }
    
};


#endif /* Dmin_h */
