//
//  Dsuite_utils.h
//  Dsuite
//
//  Created by Milan Malinsky on 02/04/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#ifndef Dsuite_utils_h
#define Dsuite_utils_h
#include <getopt.h>
#include <stdio.h>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <assert.h>
#include <time.h>
#include <regex>
#include "gzstream.h"

#define PROGRAM_BIN "Dsuite"
#define PACKAGE_BUGREPORT "milan.malinsky@unibas.ch"
#define GZIP_EXT ".gz"

#define P3isTrios2 1
#define P3isTrios1 2
#define P3isTrios0 3

using std::string;
// VCF format constant
static const int NUM_NON_GENOTYPE_COLUMNS=9;  // 8 mendatory columns + 1 column with definition of the genotype columns

double calculateOneDs(double ABBAtotal, double BABAtotal);
double* calculateThreeDs(double ABBAtotal, double BABAtotal, double BBAAtotal);
double Fd_Denom_perVariant(double p1, double p2, double p3, double pO);
double fG_Denom_perVariant(double p1, double p3a, double p3b, double pO);
double FdM_Denom_perVariant(double p1, double p2, double p3, double pO);
double normalCDF(double x);
double stringToDouble(std::string s);
std::string stripExtension(const std::string& filename);
std::vector<std::string> split2(std::string s, string delim);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<size_t> locateSet(std::vector<std::string>& sample_names, const std::vector<std::string>& set);
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in);
std::ostream* createWriter(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out);
bool file_exists(const std::string& name);
void assignTreeLevelsAndLinkToTaxa(string& treeLine, std::map<string,std::vector<int>>& taxaToLoc, std::vector<int>& levels);

// Converting numbers (int, double, size_t, and char) to string
template <typename T> std::string numToString(T i) {
    std::string ret;
    std::stringstream out;
    out << i;
    ret = out.str();
    return ret;
}

///Represents the exception for taking the median of an empty list
class median_of_empty_list_exception:public std::exception{
  virtual const char* what() const throw() {
    return "Attempt to take the median of an empty list of numbers.  "
      "The median of an empty list is undefined.";
  }
};

///Return the median of a sequence of numbers defined by the random
///access iterators begin and end.  The sequence must not be empty
///(median is undefined for an empty set).
///
///The numbers must be convertible to double.
template<class RandAccessIter> double median(RandAccessIter begin, RandAccessIter end) throw(median_of_empty_list_exception) {
  if(begin == end){ throw median_of_empty_list_exception(); }
  std::size_t size = end - begin;
  std::size_t middleIdx = size/2;
  RandAccessIter target = begin + middleIdx;
  std::nth_element(begin, target, end);

  if(size % 2 != 0){ //Odd number of elements
    return *target;
  }else{            //Even number of elements
    double a = *target;
    RandAccessIter targetNeighbor= target-1;
    std::nth_element(begin, targetNeighbor, end);
    return (a+*targetNeighbor)/2.0;
  }
}



// Print an arbitrary vector to a file
template <class T> void print_vector(T vector, std::ostream& outFile, char delim = '\t', bool endLine = true) {
    for (int i = 0; i < vector.size(); i++) {
        if (i == (vector.size()-1)) {
            if (endLine) outFile << vector[i] << std::endl;
            else outFile << vector[i];
        } else {
            outFile << vector[i] << delim;
        }
    }
}

template <class T> double vector_average(T vector) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    double average = (double)sum / (double)vector.size();
    return average;
}

template <class T> double vector_sum(T vector) {
    double sum = 0;
    for (int i = 0; i < vector.size(); i++) {
        sum += vector[i];
    }
    return sum;
}

inline void copy_except(int i, std::vector<double>& inVec, std::vector<double>& outVec) {
    std::copy(inVec.begin(), inVec.begin() + i, outVec.begin());
    std::copy(inVec.begin() + i + 1, inVec.end(), outVec.begin()+i);
    //std::cerr << "copying:" << i << " "; print_vector_stream(inVec, std::cerr);
    //std::cerr << "copied: " << i << " "; print_vector_stream(outVec, std::cerr);
}

// jackknive standard error
template <class T> double jackknive_std_err(T& vector) {
    if (vector.size() < 5) {
        throw "WARNING: Fewer than five blocks to calculate jackknife!!";
    }
    std::vector<double> jackkniveAverages;
    std::vector<double> JregionDs; JregionDs.resize(vector.size()-1);
    for (std::vector<double>::size_type i = 0; i != vector.size(); i++) {
        // std::cerr << "copying " << i << std::endl;
        copy_except(i, vector, JregionDs);
        jackkniveAverages.push_back(vector_average(JregionDs));
        JregionDs.clear(); JregionDs.resize(vector.size()-1);
    }
    double jackkniveOverallMean = vector_average(jackkniveAverages);
    double sum = 0;
    for (int i = 0; i < jackkniveAverages.size(); i++) {
        sum += std::pow((jackkniveAverages[i] - jackkniveOverallMean), 2.0);
    }
    double var = ((double)(jackkniveAverages.size()-1)/(double)jackkniveAverages.size()) * sum;
    double Dstd_err = std::sqrt(var);
    return Dstd_err;
}

class GeneralSetCounts {
public:
    GeneralSetCounts(const std::map<string, std::vector<size_t>>& setsToPosMap, const int nSamples) : overall(0) {
        for(std::map<string, std::vector<size_t>>::const_iterator it = setsToPosMap.begin(); it != setsToPosMap.end(); ++it) {
            setAltCounts[it->first] = 0; setAlleleCounts[it->first] = 0;
            setAAFs[it->first] = -1.0; setDAFs[it->first] = -1.0;
            setSizes.push_back(it->second.size());
        }
        individualsWithVariant.assign(nSamples, 0);
    };
    
    void getSetVariantCountsSimple(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    void getSetVariantCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    
    int overall;
    std::map<string,int> setAltCounts;
    std::map<string,int> setAlleleCounts; // The number of non-missing alleles for this set
    std::vector<size_t> setSizes;
    std::map<string,double> setAAFs; // Allele frequencies - alternative allele
    std::map<string,double> setDAFs; // Allele frequencies - derived allele
    std::vector<int> individualsWithVariant; // 0 homRef, 1 het, 2 homAlt
    // std::vector<int> set1individualsWithVariant; std::vector<int> set2individualsWithVariant;
    // std::vector<int> set3individualsWithVariant; std::vector<int> set4individualsWithVariant;
    
private:
    void getBasicCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
};

// Split sets for the f_G statistic
class GeneralSetCountsWithSplits : public GeneralSetCounts {
public:
    GeneralSetCountsWithSplits(const std::map<string, std::vector<size_t>>& setsToPosMap, const int nSamples) : GeneralSetCounts(setsToPosMap,nSamples) {
        for(std::map<string, std::vector<size_t>>::const_iterator it = setsToPosMap.begin(); it != setsToPosMap.end(); ++it) {
            setAAFsplit1[it->first] = -1.0; setAAFsplit2[it->first] = -1.0; setDAFsplit1[it->first] = -1.0; setDAFsplit2[it->first] = -1.0;
            setAlleleCountsSplit1[it->first] = 0; setAlleleCountsSplit2[it->first] = 0; setAltCountsSplit1[it->first] = 0; setAltCountsSplit2[it->first] = 0;
        }
    }
    std::map<string,int> setAltCountsSplit1;
    std::map<string,int> setAltCountsSplit2;
    std::map<string,double> setAAFsplit1; // Allele frequencies - alternative allele
    std::map<string,double> setAAFsplit2; //
    std::map<string,double> setDAFsplit1; // Allele frequencies - derived allele, in the complement of the set
    std::map<string,double> setDAFsplit2;
    std::map<string,int> setAlleleCountsSplit1; // The number of non-missing alleles for the complement of this set
    std::map<string,int> setAlleleCountsSplit2;
    
    void getSplitCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);

private:
    void getBasicCountsWithSplits(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
};


class TrioDinfo {
public:
    TrioDinfo() {
        ABBAtotal = 0;
        BABAtotal = 0;
        BBAAtotal = 0;
        treeArrangement = 0; BBAAarrangement = 0; DminArrangement = 0;
        regionDs.resize(3); usedVars.resize(3);
        usedVars[0] = 0; usedVars[1] = 0; usedVars[2] = 0;
        localD1num = 0; localD2num = 0; localD3num = 0;
        localD1denom = 0; localD2denom = 0; localD3denom = 0;
        F_d_denom1 = 0; F_d_denom1_reversed = 0; F_dM_denom1 = 0; F_dM_denom1_reversed = 0; F_G_denom1 = 0; F_G_denom1_reversed = 0;
        F_d_denom2 = 0; F_d_denom2_reversed = 0; F_dM_denom2 = 0; F_dM_denom2_reversed = 0; F_G_denom2 = 0; F_G_denom2_reversed = 0;
        F_d_denom3 = 0; F_d_denom3_reversed = 0; F_dM_denom3 = 0; F_dM_denom3_reversed = 0; F_G_denom3 = 0; F_G_denom3_reversed = 0;
    };
    
    // string P1; string P2; string P3;
    double ABBAtotal; double BABAtotal; double BBAAtotal;
    double D1; double D2; double D3; double D1_p; double D2_p; double D3_p;
    double F_d_denom1; double F_d_denom1_reversed; double F_dM_denom1; double F_dM_denom1_reversed; double F_G_denom1; double F_G_denom1_reversed;
    double F_d_denom2; double F_d_denom2_reversed; double F_dM_denom2; double F_dM_denom2_reversed; double F_G_denom2; double F_G_denom2_reversed;
    double F_d_denom3; double F_d_denom3_reversed; double F_dM_denom3; double F_dM_denom3_reversed; double F_G_denom3; double F_G_denom3_reversed;
    
    
    
    double localD1num; double localD2num; double localD3num;
    double localD1denom; double localD2denom; double localD3denom;
    std::vector<std::vector<double>> regionDs; // vector with three empty (double) vectors
    std::vector<int> usedVars;
    
    int treeArrangement;    // 1 - trios[i][0] and trios[i][1] are P1 and P2
                            // 2 - trios[i][0] and trios[i][2] are P1 and P2
                            // 3 - trios[i][1] and trios[i][2] are P1 and P2
    
    int BBAAarrangement;    // 1 - trios[i][0] and trios[i][1] are P1 and P2
                            // 2 - trios[i][0] and trios[i][2] are P1 and P2
                            // 3 - trios[i][1] and trios[i][2] are P1 and P2
    
    int DminArrangement;    // 1 - trios[i][0] and trios[i][1] are P1 and P2
                            // 2 - trios[i][0] and trios[i][2] are P1 and P2
                            // 3 - trios[i][1] and trios[i][2] are P1 and P2
    
    
    void assignTreeArrangement(const std::vector<int>& treeLevels, const int loc1, const int loc2, const int loc3) {
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
    
    void assignBBAAarrangement() {
        // Find which topology is in agreement with the counts of the BBAA, BABA, and ABBA patterns
        if (BBAAtotal >= BABAtotal && BBAAtotal >= ABBAtotal) {
            BBAAarrangement = P3isTrios2;
        } else if (BABAtotal >= BBAAtotal && BABAtotal >= ABBAtotal) {
            BBAAarrangement = P3isTrios1;
        } else if (ABBAtotal >= BBAAtotal && ABBAtotal >= BABAtotal) {
            BBAAarrangement = P3isTrios0;
        }
    }
    
    
    std::vector<string> makeOutVec(const std::vector<string>& trio, const bool fStats, const int arrangement) {
        
        std::vector<string> outVec; if (fStats) outVec.resize(6); else outVec.resize(5);
        
        switch (arrangement) {
                
        case P3isTrios2:
            outVec[2] = trio[2]; outVec[3] = numToString(std::fabs(D1)); outVec[4] = numToString(D1_p);
            if (D1 >= 0) {
                outVec[0] = trio[0]; outVec[1] = trio[1];
                if (fStats) {
                    double Dnum = ABBAtotal-BABAtotal;
                    outVec[5] = numToString(Dnum/F_G_denom1);
                   // outVec[6] = numToString(Dnum/F_d_denom1);
                   // outVec[7] = numToString(Dnum/F_dM_denom1);
                }
            } else {
                outVec[0] = trio[1]; outVec[1] = trio[0];
                if (fStats) {
                    double Dnum = BABAtotal-ABBAtotal;
                    outVec[5] = numToString(Dnum/F_G_denom1_reversed);
                    //outVec[6] = numToString(Dnum/F_d_denom1_reversed);
                    //outVec[7] = numToString(Dnum/F_dM_denom1_reversed);
                }
            } break;
                
        case P3isTrios1:
            outVec[2] = trio[1]; outVec[3] = numToString(std::fabs(D2)); outVec[4] = numToString(D2_p);
            if (D2 >= 0) {
                outVec[0] = trio[0]; outVec[1] = trio[2];
                if (fStats) {
                    double Dnum = ABBAtotal - BBAAtotal;
                    outVec[5] = numToString(Dnum/F_G_denom2);
                   // outVec[6] = numToString(Dnum/F_d_denom2);
                   // outVec[7] = numToString(Dnum/F_dM_denom2);
                }
            } else {
                outVec[0] = trio[2]; outVec[1] = trio[0];
                if (fStats) {
                    double Dnum = BBAAtotal - ABBAtotal;
                    outVec[5] = numToString(Dnum/F_G_denom2_reversed);
                   // outVec[6] = numToString(Dnum/F_d_denom2_reversed);
                   // outVec[7] = numToString(Dnum/F_dM_denom2_reversed);
                }
            } break;
                
        case P3isTrios0:
            outVec[2] = trio[0]; outVec[3] = numToString(std::fabs(D3)); outVec[4] = numToString(D3_p);
            if (D3 >= 0) {
                outVec[0] = trio[2]; outVec[1] = trio[1];
                if (fStats) {
                    double Dnum = BBAAtotal - BABAtotal;
                    outVec[5] = numToString(Dnum/F_G_denom3);
                   // outVec[6] = numToString(Dnum/F_d_denom3);
                   // outVec[7] = numToString(Dnum/F_dM_denom3);
                }
            } else {
                outVec[0] = trio[1]; outVec[1] = trio[2];
                if (fStats) {
                    double Dnum = BABAtotal - BBAAtotal;
                    outVec[5] = numToString(Dnum/F_G_denom3_reversed);
                   // outVec[6] = numToString(Dnum/F_d_denom3_reversed);
                   // outVec[7] = numToString(Dnum/F_dM_denom3_reversed);
                }
            } break;
                
        }
        return outVec;
    }
    
    void assignDminArrangement() {
        
        if (std::fabs(D1) <= std::fabs(D2) && std::fabs(D1) <= std::fabs(D3)) {
            DminArrangement = P3isTrios2;
        } else if (std::fabs(D2) <= std::fabs(D1) && std::fabs(D2) <= std::fabs(D3)) { // (P3 == S2)
            DminArrangement = P3isTrios1;
        } else if (std::fabs(D3) <= std::fabs(D1) && std::fabs(D3) <= std::fabs(D2)) { // (P3 == S1)
            DminArrangement = P3isTrios0;
        }
        
    }
    
    
    
    void calculateFinalDs() throw(const char*) {
        double* Ds = calculateThreeDs(ABBAtotal, BABAtotal, BBAAtotal);
        D1 = Ds[0]; D2 = Ds[1]; D3 = Ds[2];
        
        // Get the standard error values:
        double D1stdErr = jackknive_std_err(regionDs[0]); double D2stdErr = jackknive_std_err(regionDs[1]);
        double D3stdErr = jackknive_std_err(regionDs[2]);
        //std::cerr << "Here: " << regionDs[2][0] << std::endl;
        // Get the Z-scores
        double D1_Z = std::fabs(D1)/D1stdErr; double D2_Z = std::fabs(D2)/D2stdErr;
        double D3_Z = std::fabs(D3)/D3stdErr;
        // And p-values
        D1_p = 1 - normalCDF(D1_Z); D2_p = 1 - normalCDF(D2_Z);
        D3_p = 1 - normalCDF(D3_Z);
    }
    
    void addRegionDs(const int arrangement) {
        switch (arrangement) {
            case P3isTrios2: regionDs[0].push_back(localD1num/localD1denom); localD1num = 0; localD1denom = 0; usedVars[0] = 0; break;
            case P3isTrios1: regionDs[1].push_back(localD2num/localD2denom); localD2num = 0; localD2denom = 0; usedVars[1] = 0; break;
            case P3isTrios0: regionDs[2].push_back(localD3num/localD3denom); localD3num = 0; localD3denom = 0; usedVars[2] = 0; break;
        }
    }
    
};


#endif /* Dsuite_utils_h */

