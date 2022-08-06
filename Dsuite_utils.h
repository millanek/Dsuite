//
//  Dsuite_utils.h
//  Dsuite
//
//  Created by Milan Malinsky on 02/04/2019.
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
#include <algorithm>
#include <limits>
#include <random>
#include <list>
#include <cstdint>
#include "gzstream.h"
#include "kstest.h"

#define PROGRAM_BIN "Dsuite"
#define PACKAGE_BUGREPORT "milan.malinsky@iee.unibe.ch"
#define GZIP_EXT ".gz"
#define F4HEADER "f4-ratio"
#define ploidy 2

#define CHECK_TREE_ERROR_MSG "It seems that this species is in the SETS.txt file but can't be found in the tree. Please check the spelling and completeness of your tree file."

#define P3isTrios2_Dpositive 1      // 1 - trios[i][0] and trios[i][1] are P1 and P2; D >= 0
#define P3isTrios2_Dnegative 2      // 2 - trios[i][0] and trios[i][1] are P1 and P2; D < 0
#define P3isTrios1_Dpositive 3      // 3 - trios[i][0] and trios[i][2] are P1 and P2; D >= 0
#define P3isTrios1_Dnegative 4      // 4 - trios[i][0] and trios[i][2] are P1 and P2; D < 0
#define P3isTrios0_Dpositive 5      // 5 - trios[i][2] and trios[i][1] are P1 and P2; D >= 0
#define P3isTrios0_Dnegative 6      // 6 - trios[i][2] and trios[i][1] are P1 and P2; D < 0

#define P3isTrios2 7    // 7 - trios[i][0] and trios[i][1] are P1 and P2;
#define P3isTrios1 8    // 8 - trios[i][0] and trios[i][2] are P1 and P2;
#define P3isTrios0 9    // 9 - trios[i][1] and trios[i][2] are P1 and P2;

#define ABBAvector 0
#define BABAvector 1
#define BBAAvector 2


#define OutgroupNotRequired 0
#define OutgroupRequired 1

#define LikelihoodsProbabilitiesAbsent 0
#define LikelihoodsProbabilitiesGP 1
#define LikelihoodsProbabilitiesGL 2
#define LikelihoodsProbabilitiesPL 3

#define AncestralAlleleMissing -1
#define AncestralAlleleRef 0
#define AncestralAlleleAlt 1

using std::string;
// VCF format constant
static const int NUM_NON_GENOTYPE_COLUMNS=9;  // 8 mendatory columns + 1 column with definition of the genotype columns

double calculateOneDs(double ABBAtotal, double BABAtotal);
double* calculateThreeDs(double ABBAtotal, double BABAtotal, double BBAAtotal);
double Fd_Denom_perVariant(double p1, double p2, double p3, double pO);
double fG_Denom_perVariant(double p1, double p3a, double p3b, double pO);
double FdM_Denom_perVariant(double p1, double p2, double p3, double pO);
long double normalCDF(double x);
double stringToDouble(std::string s);
std::string stripExtension(const std::string& filename);
std::vector<std::string> split2(std::string s, string delim);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<double> splitToDouble(const std::string &s, char delim);
std::vector<size_t> locateSet(std::vector<std::string>& sample_names, const std::vector<std::string>& set);
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in);
std::ostream* createWriter(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out);
bool file_exists(const std::string& name);
void assignSplits01FromAlleleFrequency(const double p, double& splitA, double& splitB);

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
template<class RandAccessIter> double median(RandAccessIter begin, RandAccessIter end) {
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

inline unsigned nChoosek( unsigned n, unsigned k )
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;
    
    int result = n;
    for( int i = 2; i <= k; ++i ) {
        result *= (n-i+1);
        result /= i;
    }
    return result;
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
    GeneralSetCounts(const std::map<string, std::vector<size_t>>& setsToPosMap, const int nSamples) : overall(0), averageAAF(-1.0), averageDAF(-1.0),  likelihoodsProbabilitiesType(LikelihoodsProbabilitiesAbsent), AAint(AncestralAlleleMissing) {
        for(std::map<string, std::vector<size_t>>::const_iterator it = setsToPosMap.begin(); it != setsToPosMap.end(); ++it) {
            setAltCounts[it->first] = 0; setAlleleCounts[it->first] = 0; setAlleleProbCounts[it->first] = 0;
            setAAFs[it->first] = -1.0; setDAFs[it->first] = -1.0;
            setAAFsFromLikelihoods[it->first] = -1.0; setDAFsFromLikelihoods[it->first] = -1.0;
            setPoolAAFs[it->first] = -1.0; setPoolDAFs[it->first] = -1.0;
            setSizes.push_back(it->second.size()); setCorrectionFactors[it->first] = -1.0;
            setHWEpriorsFromAAFfromGT[it->first].assign(3, -1.0);
            setHWEpriorsFromDAFfromGT[it->first].assign(3, -1.0);
            std::vector<int> thisSetGenotypes; setGenotypes[it->first] = thisSetGenotypes;
            std::vector<double> thisSetIndividualAFs; setIndividualAFs[it->first] = thisSetIndividualAFs;
        }
        individualsWithVariant.assign(nSamples, 0);
        individualPoolAAFs.assign(nSamples, -1.0);
    };
    
    void getSetVariantCountsSimple(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    void getSetVariantCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    
    int checkForGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& vcfLineFields);
    int findADtagPosition(const std::vector<std::string>& vcfLineFields);
    
    void getAFsFromGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& genotypeFields, const std::map<size_t, string>& posToSpeciesMap, const int likelihoodsOrProbabilitiesTagPosition);
    void getAFsFromADtag(const std::vector<std::string>& genotypeFields, const std::map<string, std::vector<size_t>>& setsToPosMap, const int ADTagPosition, const int minDepth);
    
    int overall; int AAint;
    std::map<string,std::vector<int>> setGenotypes;
    std::map<string,std::vector<double>> setIndividualAFs;
    std::map<string,int> setAltCounts;
    std::map<string,int> setAlleleCounts; // The number of non-missing alleles for this set
    std::map<string,int> setAlleleProbCounts; // The number of non-missing alleles for this set in terms of likelihoods/probabilities
    std::vector<double> individualPoolAAFs;  // Allele frequency for each individual pool estimated from Allelic Depth (AD tag in VCF) - for pool-seq data
    std::map<string,double> setPoolAAFs; // The above individual pool values are then averaged if multiple pools form a set (i.e., a population or species)
    std::map<string,double> setPoolDAFs;
    
    std::vector<size_t> setSizes; std::map<string,double> setCorrectionFactors;
    std::map<string,double> setAAFs; double averageAAF;     // Allele frequencies - alternative allele
    std::map<string,double> setDAFs; double averageDAF;     // Allele frequencies - derived allele
    std::map<string,double> setAAFsFromLikelihoods; double averageAAFFromLikelihoods;   // Allele frequencies - alternative allele
    std::map<string,double> setDAFsFromLikelihoods; double averageDAFFromLikelihoods;   // Allele frequencies - derived allele
    std::vector<int> individualsWithVariant; // 0 homRef, 1 het, 2 homAlt
    int likelihoodsProbabilitiesType;
    // std::vector<int> set1individualsWithVariant; std::vector<int> set2individualsWithVariant;
    // std::vector<int> set3individualsWithVariant; std::vector<int> set4individualsWithVariant;
    

    int returnFormatTagPosition(std::vector<std::string>& format, const std::string& tag);
    void setHWEpriorsFromAFfromGT();
    std::vector<double> probabilitiesFromLikelihoods(const std::vector<double>& thisLikelihoods, const string& species);
    std::map<string,std::vector<double> > setHWEpriorsFromAAFfromGT;
    std::map<string,std::vector<double> > setHWEpriorsFromDAFfromGT;
    
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
            
            setAAFsplit1fromLikelihoods[it->first] = -1.0; setAAFsplit2fromLikelihoods[it->first] = -1.0; setDAFsplit1fromLikelihoods[it->first] = -1.0;
            setDAFsplit2fromLikelihoods[it->first] = -1.0; setAlleleCountsSplit1fromLikelihoods[it->first] = 0;
            setAlleleCountsSplit2fromLikelihoods[it->first] = 0;
            
            setPoolAAFsplit1[it->first] = -1.0; setPoolAAFsplit2[it->first] = -1.0;
            setPoolDAFsplit1[it->first] = -1.0; setPoolDAFsplit2[it->first] = -1.0;

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
    
    

    std::map<string,double> setAAFsplit1fromLikelihoods; // Allele frequencies - alternative allele
    std::map<string,double> setAAFsplit2fromLikelihoods; //
    std::map<string,double> setDAFsplit1fromLikelihoods; // Allele frequencies - derived allele, in the complement of the set
    std::map<string,double> setDAFsplit2fromLikelihoods;
    std::map<string,int> setAlleleCountsSplit1fromLikelihoods; // The number of non-missing alleles for the complement of this set
    std::map<string,int> setAlleleCountsSplit2fromLikelihoods;
    
    std::map<string,double> setPoolAAFsplit1;
    std::map<string,double> setPoolDAFsplit1;
    std::map<string,double> setPoolAAFsplit2; 
    std::map<string,double> setPoolDAFsplit2;
    
    
    
    void getSplitCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    void getSplitCountsNew(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    
    void getAFsFromADtagWithSplits(const std::vector<std::string>& genotypeFields, const std::map<string, std::vector<size_t>>& setsToPosMap, const int ADTagPosition, const int minDepth);
    void getAFsFromGenotypeLikelihoodsOrProbabilitiesWithSplits(const std::vector<std::string>& genotypeFields, const std::map<size_t, string>& posToSpeciesMap, const int likelihoodsOrProbabilitiesTagPosition);

private:
    void getBasicCountsWithSplits(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    void getBasicCountsWithSplitsNew(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
};

class TrioDinfo {
public:
    TrioDinfo() {
        ABBAtotal = 0;
        BABAtotal = 0;
        BBAAtotal = 0;
        treeArrangement = 0; BBAAarrangement = 0; DminArrangement = 0;
        regionDs.resize(3); usedVars.resize(3); totalUsedVars.resize(3);
        linearStrongABBApos.resize(3);  linearStrongBABApos.resize(3);
        usedVars[0] = 0; usedVars[1] = 0; usedVars[2] = 0;
        localD1num = 0; localD2num = 0; localD3num = 0;
        localD1denom = 0; localD2denom = 0; localD3denom = 0;
        F_d_denom1 = 0; F_d_denom1_reversed = 0; F_dM_denom1 = 0; F_dM_denom1_reversed = 0; F_G_denom1 = 0; F_G_denom1_reversed = 0;
        F_d_denom2 = 0; F_d_denom2_reversed = 0; F_dM_denom2 = 0; F_dM_denom2_reversed = 0; F_G_denom2 = 0; F_G_denom2_reversed = 0;
        F_d_denom3 = 0; F_d_denom3_reversed = 0; F_dM_denom3 = 0; F_dM_denom3_reversed = 0; F_G_denom3 = 0; F_G_denom3_reversed = 0;
    };
    
    // string P1; string P2; string P3;
    double ABBAtotal; double BABAtotal; double BBAAtotal;
    double D1; double D2; double D3; long double D1_p; long double D2_p; long double D3_p; double D1_Z; double D2_Z; double D3_Z;
    double F_d_denom1; double F_d_denom1_reversed; double F_dM_denom1; double F_dM_denom1_reversed; double F_G_denom1; double F_G_denom1_reversed;
    double F_d_denom2; double F_d_denom2_reversed; double F_dM_denom2; double F_dM_denom2_reversed; double F_G_denom2; double F_G_denom2_reversed;
    double F_d_denom3; double F_d_denom3_reversed; double F_dM_denom3; double F_dM_denom3_reversed; double F_G_denom3; double F_G_denom3_reversed;
    
    std::vector<std::vector<int>> linearStrongABBApos; // positions of strong (> 0.5) ABBA for the three tree orientations
    std::vector<std::vector<int>> linearStrongBABApos; // positions of strong (> 0.5) BABA for the three tree orientations
    
    double KSpvalForStrongSites;
    double KSpvalForStrongSitesOneSample;
    
    
    
    double localD1num; double localD2num; double localD3num;
    double localD1denom; double localD2denom; double localD3denom;
    std::vector<std::vector<double>> regionDs; // vector with three empty (double) vectors
    std::vector<int> usedVars; // Used vars for local windows
    std::vector<int> totalUsedVars; // Used vars for local windows
    
    
    // 1 - trios[i][0] and trios[i][1] are P1 and P2; D >= 0
    // 2 - trios[i][0] and trios[i][1] are P1 and P2; D < 0
    // 3 - trios[i][0] and trios[i][2] are P1 and P2; D >= 0
    // 4 - trios[i][0] and trios[i][2] are P1 and P2; D < 0
    // 5 - trios[i][1] and trios[i][2] are P1 and P2; D >= 0
    // 6 - trios[i][1] and trios[i][2] are P1 and P2; D < 0
    int treeArrangement; int BBAAarrangement; int DminArrangement;    
    
    
    int assignTreeArrangement(const std::vector<int>& treeLevels, const int loc1, const int loc2, const int loc3) {
        int midLoc = std::max(std::min(loc1,loc2), std::min(std::max(loc1,loc2),loc3));
        int arrangement = 0;
        if (midLoc == loc1) {
            if (loc2 < loc1) {
                int m1 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc1);
                int m2 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc3);
                if (m1 < m2) { if (D2 >= 0) { arrangement = P3isTrios1_Dpositive; } else { arrangement = P3isTrios1_Dnegative;}}
                else { if (D1 >= 0) { arrangement = P3isTrios2_Dpositive; } else { arrangement = P3isTrios2_Dnegative; } }
            } else {
                int m1 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc1);
                int m2 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc2);
                if (m1 < m2) { if (D1 >= 0) { arrangement = P3isTrios2_Dpositive; } else { arrangement = P3isTrios2_Dnegative; }}
                else { if (D2 >= 0) arrangement = P3isTrios1_Dpositive; else arrangement = P3isTrios1_Dnegative;}
            }
        } else if (midLoc == loc2) {
            if (loc1 < loc2) {
                int m1 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc2);
                int m2 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc3);
                if (m1 < m2) { if (D3 >= 0) arrangement = P3isTrios0_Dpositive; else arrangement = P3isTrios0_Dnegative; }
                else { if (D1 >= 0) { arrangement = P3isTrios2_Dpositive; } else { arrangement = P3isTrios2_Dnegative; }}
            } else {
                int m1 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc2);
                int m2 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc1);
                if (m1 < m2) { if (D1 >= 0) { arrangement = P3isTrios2_Dpositive;} else {arrangement = P3isTrios2_Dnegative;} }
                else { if (D3 >= 0) arrangement = P3isTrios0_Dpositive; else arrangement = P3isTrios0_Dnegative; }
            }
        } else if (midLoc == loc3) {
            if (loc1 < loc3) {
                int m1 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc3);
                int m2 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc2);
                if (m1 < m2) { if (D3 >= 0) arrangement = P3isTrios0_Dpositive; else arrangement = P3isTrios0_Dnegative; }
                else { if (D2 >= 0) arrangement = P3isTrios1_Dpositive; else arrangement = P3isTrios1_Dnegative;}
                
            } else {
                int m1 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc3);
                int m2 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc1);
                if (m1 < m2) { if (D2 >= 0) arrangement = P3isTrios1_Dpositive; else arrangement = P3isTrios1_Dnegative;}
                else { if (D3 >= 0) arrangement = P3isTrios0_Dpositive; else arrangement = P3isTrios0_Dnegative; }
            }
        }
        return arrangement;
    }
    
    void assignBBAAarrangement() {
        // Find which topology is in agreement with the counts of the BBAA, BABA, and ABBA patterns
        if (BBAAtotal >= BABAtotal && BBAAtotal >= ABBAtotal) {
            if (D1 >= 0) BBAAarrangement = P3isTrios2_Dpositive; else BBAAarrangement = P3isTrios2_Dnegative;
        } else if (BABAtotal >= BBAAtotal && BABAtotal >= ABBAtotal) {
            if (D2 >= 0) BBAAarrangement = P3isTrios1_Dpositive; else BBAAarrangement = P3isTrios1_Dnegative;
        } else if (ABBAtotal >= BBAAtotal && ABBAtotal >= BABAtotal) {
            if (D3 >= 0) BBAAarrangement = P3isTrios0_Dpositive; else BBAAarrangement = P3isTrios0_Dnegative;
        }
    }
    
    
    std::vector<string> makeOutVec(const std::vector<string>& trio, const bool fStats, const bool KS_test, const int arrangement) {
        
        int numCols = 9; if (fStats) numCols++; if(KS_test) numCols++;
        std::vector<string> outVec; outVec.resize(numCols);
        int patternsI = 6; if (fStats) patternsI++; if(KS_test) patternsI++; // Where will the BBAA, ABBA, etc. counts be put
        
        
        if (std::fabs(D1_p) < 2.3e-16) { D1_p = 2.3e-16; }
        if (std::fabs(D2_p) < 2.3e-16) { D2_p = 2.3e-16; }
        if (std::fabs(D3_p) < 2.3e-16) { D3_p = 2.3e-16; }
        
        double f4ratio; double Dnum;
        
        switch (arrangement) {
                
        case P3isTrios2_Dpositive:
            outVec[0] = trio[0]; outVec[1] = trio[1]; outVec[2] = trio[2];
            outVec[3] = numToString(D1); outVec[4] = numToString(D1_Z);
            outVec[5] = numToString(D1_p); outVec[patternsI] = numToString(BBAAtotal);
            outVec[patternsI+1] = numToString(ABBAtotal); outVec[patternsI+2] = numToString(BABAtotal);
            if (fStats) { Dnum = ABBAtotal-BABAtotal; f4ratio = Dnum/F_G_denom1; }
            break;
                
        case P3isTrios2_Dnegative:
            outVec[0] = trio[1]; outVec[1] = trio[0]; outVec[2] = trio[2];
            outVec[3] = numToString(std::fabs(D1)); outVec[4] = numToString(D1_Z);
            outVec[5] = numToString(D1_p); outVec[patternsI] = numToString(BBAAtotal);
            outVec[patternsI+1] = numToString(BABAtotal); outVec[patternsI+2] = numToString(ABBAtotal);
            if (fStats) { Dnum = BABAtotal-ABBAtotal; f4ratio = Dnum/F_G_denom1_reversed; }
            break;
                
        case P3isTrios1_Dpositive:
            outVec[0] = trio[0]; outVec[1] = trio[2];
            outVec[2] = trio[1]; outVec[3] = numToString(std::fabs(D2)); outVec[4] = numToString(D2_Z);
            outVec[5] = numToString(D2_p); outVec[patternsI] = numToString(BABAtotal);
            outVec[patternsI+1] = numToString(ABBAtotal); outVec[patternsI+2] = numToString(BBAAtotal);
            if (fStats) { Dnum = ABBAtotal - BBAAtotal; f4ratio = Dnum/F_G_denom2; }
            break;
        
        case P3isTrios1_Dnegative:
            outVec[0] = trio[2]; outVec[1] = trio[0];
            outVec[2] = trio[1]; outVec[3] = numToString(std::fabs(D2)); outVec[4] = numToString(D2_Z);
            outVec[5] = numToString(D2_p); outVec[patternsI] = numToString(BABAtotal);
            outVec[patternsI+1] = numToString(BBAAtotal); outVec[patternsI+2] = numToString(ABBAtotal);
            if (fStats) { Dnum = BBAAtotal - ABBAtotal; f4ratio = Dnum/F_G_denom2_reversed; }
            break;
                
        case P3isTrios0_Dpositive:
            outVec[0] = trio[2]; outVec[1] = trio[1];
            outVec[2] = trio[0]; outVec[3] = numToString(std::fabs(D3)); outVec[4] = numToString(D3_Z);
            outVec[5] = numToString(D3_p); outVec[patternsI] = numToString(ABBAtotal);
            outVec[patternsI+1] = numToString(BBAAtotal); outVec[patternsI+2] = numToString(BABAtotal);
            if (fStats) { Dnum = BBAAtotal - BABAtotal; f4ratio = Dnum/F_G_denom3; }
            break;
        
        case P3isTrios0_Dnegative:
            outVec[0] = trio[1]; outVec[1] = trio[2];
            outVec[2] = trio[0]; outVec[3] = numToString(std::fabs(D3)); outVec[4] = numToString(D3_Z);
            outVec[5] = numToString(D3_p); outVec[patternsI] = numToString(ABBAtotal);
            outVec[patternsI+1] = numToString(BABAtotal); outVec[patternsI+2] = numToString(BBAAtotal);
            if (fStats) { Dnum = BABAtotal - BBAAtotal; f4ratio = Dnum/F_G_denom3_reversed; }
            break;
                
        }
        
        /* investigating rare cases of unexpected f4-ratio values
        if (f4ratio < 0 || f4ratio > 1) {
            std::cerr << "trio[0]: " << trio[0] << "trio[1]: " << trio[1] << "trio[2]: " << trio[2] << std::endl;
            std::cerr << "Arrangement: " << arrangement << std::endl;
            std::cerr << "D1: " << D1 << " D2: " << D2  << " D3: " << D3 << std::endl;
            std::cerr << "ABBAtotal: " << ABBAtotal << "BBAAtotal: " << BBAAtotal  << "BABAtotal: " << BABAtotal << std::endl;
            std::cerr << "ABBAtotal-BABAtotal: " << ABBAtotal-BABAtotal << "ABBAtotal-BBAAtotal: " << ABBAtotal- BBAAtotal  << "BBAAtotal-BABAtotal: " << BBAAtotal-BABAtotal << std::endl;
            std::cerr << "F_G_denom1: " << F_G_denom1 << "; F_G_denom1_reversed: " << F_G_denom1_reversed << std::endl;
            std::cerr << "F_G_denom2: " << F_G_denom2 << "; F_G_denom2_reversed: " << F_G_denom2_reversed << std::endl;
            std::cerr << "F_G_denom3: " << F_G_denom3 << "; F_G_denom3_reversed: " << F_G_denom3_reversed << std::endl;
            std::cerr << "f4ratio: " << f4ratio << std::endl;
        }
        */
        
        
        if (fStats) { // For now just bounding those unexpected cases to the 0 to 1 range
            if (f4ratio < 0) { f4ratio = 0; }
            if (f4ratio > 1) { f4ratio = 1; }
            outVec[6] = numToString(f4ratio);
        }
        if (KS_test) {
            if (fStats) outVec[7] = numToString(KSpvalForStrongSitesOneSample);
            else outVec[6] = numToString(KSpvalForStrongSitesOneSample);
        }
        
        return outVec;
    }
    
    void assignDminArrangement() {
        
        if (std::fabs(D1) <= std::fabs(D2) && std::fabs(D1) <= std::fabs(D3)) {
            if (D1 >= 0) DminArrangement = P3isTrios2_Dpositive; else DminArrangement = P3isTrios2_Dnegative;
        } else if (std::fabs(D2) <= std::fabs(D1) && std::fabs(D2) <= std::fabs(D3)) { // (P3 == S2)
            if (D2 >= 0) DminArrangement = P3isTrios1_Dpositive; else DminArrangement = P3isTrios1_Dnegative;
        } else if (std::fabs(D3) <= std::fabs(D1) && std::fabs(D3) <= std::fabs(D2)) { // (P3 == S1)
            if (D3 >= 0) DminArrangement = P3isTrios0_Dpositive; else DminArrangement = P3isTrios0_Dnegative;
        }
        
    }
    
    
    
    void calculateFinalDs() {
        double* Ds = calculateThreeDs(ABBAtotal, BABAtotal, BBAAtotal);
        D1 = Ds[0]; D2 = Ds[1]; D3 = Ds[2];
        
        // Get the standard error values:
        double D1stdErr = jackknive_std_err(regionDs[0]); double D2stdErr = jackknive_std_err(regionDs[1]);
        double D3stdErr = jackknive_std_err(regionDs[2]);
        //std::cerr << "Here: " << regionDs[2][0] << std::endl;
        // Get the Z-scores
        D1_Z = std::fabs(D1)/D1stdErr; D2_Z = std::fabs(D2)/D2stdErr; D3_Z = std::fabs(D3)/D3stdErr;
        // And p-values
        D1_p = 2 * (1 - normalCDF(D1_Z)); D2_p = 2 * (1 - normalCDF(D2_Z));
        D3_p = 2 * (1 - normalCDF(D3_Z));
    }
    
    void testIfStrongSitesUniformlyDistributed() {
        std::random_device rd;     // only used once to initialise (seed) engine
        std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
        
      // std::cerr << "Here: .. treeArrangement: " << BBAAarrangement << std::endl;
        
        std::list<int64_t> linearStrongPosList;
        std::vector<int> linearStrongPosVector;
        
        int usedVarsThisArrangment = 0;
        switch (BBAAarrangement) {
            case P3isTrios2_Dpositive:
                usedVarsThisArrangment = totalUsedVars[0];
                //linearStrongPosList.assign(linearStrongABBApos[0].begin(),linearStrongABBApos[0].end());
                linearStrongPosVector.assign(linearStrongABBApos[0].begin(),linearStrongABBApos[0].end());
                break;
            case P3isTrios2_Dnegative:
                usedVarsThisArrangment = totalUsedVars[0];
                //linearStrongPosList.assign(linearStrongBABApos[0].begin(),linearStrongBABApos[0].end());
                linearStrongPosVector.assign(linearStrongBABApos[0].begin(),linearStrongBABApos[0].end());
                break;
            case P3isTrios1_Dpositive:
                usedVarsThisArrangment = totalUsedVars[1];
                //linearStrongPosList.assign(linearStrongABBApos[1].begin(),linearStrongABBApos[1].end());
                linearStrongPosVector.assign(linearStrongABBApos[1].begin(),linearStrongABBApos[1].end());
                break;
            case P3isTrios1_Dnegative:
                usedVarsThisArrangment = totalUsedVars[1];
                //linearStrongPosList.assign(linearStrongBABApos[1].begin(),linearStrongBABApos[1].end());
                linearStrongPosVector.assign(linearStrongBABApos[1].begin(),linearStrongBABApos[1].end());
                break;
            case P3isTrios0_Dpositive:
                usedVarsThisArrangment = totalUsedVars[2];
                //linearStrongPosList.assign(linearStrongABBApos[2].begin(),linearStrongABBApos[2].end());
                linearStrongPosVector.assign(linearStrongABBApos[2].begin(),linearStrongABBApos[2].end());
                break;
            case P3isTrios0_Dnegative:
                usedVarsThisArrangment = totalUsedVars[2];
                //linearStrongPosList.assign(linearStrongBABApos[2].begin(),linearStrongBABApos[2].end());
                linearStrongPosVector.assign(linearStrongBABApos[2].begin(),linearStrongBABApos[2].end());
                break;
        }
        
        if(linearStrongPosVector.size() < 2) {
           // KSpvalForStrongSites = 1;
            KSpvalForStrongSitesOneSample = 1;
        } else {
        
         //   std::uniform_int_distribution<int> uniABBA(0,usedVarsThisArrangment); // guaranteed unbiased
         //   int numUniformSamples = (int)linearStrongPosList.size(); if (numUniformSamples < 10000) { numUniformSamples = usedVarsThisArrangment/4; }
         //   std::list<int64_t> uniABBAvals; std::vector<int64_t> uniABBAvalsVec; for (int i = 0; i < numUniformSamples; i++) { int64_t r = uniABBA(rng); uniABBAvals.push_back(r); uniABBAvalsVec.push_back(r); }
            
            //std::cerr << "Here: .. uniABBAvals.size(): " << uniABBAvals.size() << std::endl;
           // std::cerr << "Here: .. usedVarsThisArrangment: " << usedVarsThisArrangment << std::endl;
           // std::cerr << "Here: .. linearStrongPosList.size(): " << linearStrongPosList.size() << std::endl;
            
            
            
           // std::cerr << "Here: .. linearStrongPosList.size(): " << linearStrongPosList.size() << std::endl;
            
            
          //  KSpvalForStrongSites = ks_test(uniABBAvals, linearStrongPosList, std::cerr, false);
        
            
           // std::vector<int> linearStrongPosListVec(linearStrongPosList.begin(), linearStrongPosList.end());
            std::vector<double> linearStrongPosVector0to1(linearStrongPosVector.size(),0.0);
            for (int i = 0; i < linearStrongPosVector0to1.size(); i++) {
                linearStrongPosVector0to1[i] = (double)linearStrongPosVector[i]/(double)usedVarsThisArrangment;
            }
            
            KSpvalForStrongSitesOneSample = ks_test_of_uniformity(linearStrongPosVector0to1, std::cerr, false);
            
          /*  std::cerr << "KSpvalForStrongSites: " << KSpvalForStrongSites << std::endl;
            std::cerr << "KSpvalForStrongSitesOneSample: " << KSpvalForStrongSitesOneSample << std::endl;
            std::cerr << std::endl;
            
            //if (KSpvalForStrongSites > 0.99 && KSpvalForStrongSitesOneSample < 0.00001) {
            if(usedVarsThisArrangment < 10 || linearStrongPosList.size() < 5) {
                std::cerr << "KSpvalForStrongSites: " << KSpvalForStrongSites << std::endl;
                std::cerr << "KSpvalForStrongSitesOneSample: " << KSpvalForStrongSitesOneSample << std::endl;
                std::cerr << "usedVarsThisArrangment: " << usedVarsThisArrangment << std::endl;
                std::cerr << "uniABBAvals: " << std::endl;
                print_vector(uniABBAvalsVec, std::cerr, ',');
                std::cerr << "linearStrongPosListVec: " << std::endl;
                print_vector(linearStrongPosListVec, std::cerr, ',');
                std::cerr << "linearStrongPosVector0to1: " << std::endl;
                print_vector(linearStrongPosVector0to1, std::cerr, ',');
                exit(1);
            }
            */
           // std::cerr << "Here: .. KSpvalForStrongSites: " << KSpvalForStrongSites << std::endl;
            //double BABApval = ks_test(uniBABAvals, linearBABApos, std::cerr);
       }
    }
    
    
    
    void addRegionDs(const int arrangement) {
        switch (arrangement) {
            case P3isTrios2: if(localD1denom > 0) regionDs[0].push_back(localD1num/localD1denom); localD1num = 0; localD1denom = 0; usedVars[0] = 0; break;
            case P3isTrios1: if(localD2denom > 0) regionDs[1].push_back(localD2num/localD2denom); localD2num = 0; localD2denom = 0; usedVars[1] = 0; break;
            case P3isTrios0: if(localD3denom > 0) regionDs[2].push_back(localD3num/localD3denom); localD3num = 0; localD3denom = 0; usedVars[2] = 0; break;
        }
    }
    
};

//  P3 \     / P2
//      -----
//  P4 /     \ P1


/*
 Quartets are unrooted, so there are three arrangements
 (P1,P2)(P3,P4)   // BBAA (AABB) highest
 (P1,P3)(P2,P4)   // BABA (ABAB) highest
 (P1,P4)(P2,P3)   // ABBA (BAAB) highest
 
 */
class QuartetDinfo: public TrioDinfo {
public:
    
    int assignQuartetTreeArrangement(const std::vector<int>& treeLevels, const int loc1, const int loc2, const int loc3, const int loc4) {
        int firstThreeArranged = assignTreeArrangement(treeLevels, loc1, loc2, loc3);
        int withFourthArranged = 0; int overallTreeArrangment = 0;
        switch (firstThreeArranged) {
            case P3isTrios2:
                withFourthArranged = assignTreeArrangement(treeLevels, loc1, loc2, loc4);
                switch (withFourthArranged) {
                    case P3isTrios2: overallTreeArrangment = P3isTrios2; break;
                    case P3isTrios1: overallTreeArrangment = P3isTrios0; break;
                    case P3isTrios0: overallTreeArrangment = P3isTrios1; break;
                } break;
            case P3isTrios1:
                withFourthArranged = assignTreeArrangement(treeLevels, loc1, loc3, loc4);
                switch (withFourthArranged) {
                    case P3isTrios2: overallTreeArrangment = P3isTrios1; break;
                    case P3isTrios1: overallTreeArrangment = P3isTrios0; break;
                    case P3isTrios0: overallTreeArrangment = P3isTrios2; break;
                } break;
            case P3isTrios0:
            withFourthArranged = assignTreeArrangement(treeLevels, loc2, loc3, loc4);
            switch (withFourthArranged) {
                case P3isTrios2: overallTreeArrangment = P3isTrios0; break;
                case P3isTrios1: overallTreeArrangment = P3isTrios1; break;
                case P3isTrios0: overallTreeArrangment = P3isTrios2; break;
            } break;
        }
        return overallTreeArrangment;
    }
    
    
    std::vector<string> makeOutVec(const std::vector<string>& quartet, const bool fStats, const int arrangement) {
        
        std::vector<string> outVec; if (fStats) outVec.resize(11); else outVec.resize(10);
        int patternsI; if (fStats) patternsI = 8; else patternsI = 7; // Where will be put the BBAA, ABBA, etc. counts
        
        switch (arrangement) {
                
        case P3isTrios2:
            outVec[2] = quartet[2]; outVec[3] = quartet[3];
            outVec[4] = numToString(std::fabs(D1)); outVec[5] = numToString(D1_Z);
            outVec[6] = numToString(D1_p); outVec[patternsI] = numToString(BBAAtotal);
            if (D1 >= 0) {
                outVec[0] = quartet[0]; outVec[1] = quartet[1];
                outVec[patternsI+1] = numToString(ABBAtotal); outVec[patternsI+2] = numToString(BABAtotal);
                if (fStats) {
                    double Dnum = ABBAtotal-BABAtotal;
                    outVec[7] = numToString(Dnum/F_G_denom1);
                   
                }
            } else {
                outVec[0] = quartet[1]; outVec[1] = quartet[0];
                outVec[patternsI+1] = numToString(BABAtotal); outVec[patternsI+2] = numToString(ABBAtotal);
                if (fStats) {
                    double Dnum = BABAtotal-ABBAtotal;
                    outVec[7] = numToString(Dnum/F_G_denom1_reversed);
                }
            } break;
                
        case P3isTrios1:
            outVec[2] = quartet[1]; outVec[3] = quartet[3];
            outVec[4] = numToString(std::fabs(D2)); outVec[5] = numToString(D2_Z);
            outVec[6] = numToString(D2_p); outVec[patternsI] = numToString(BABAtotal);
            if (D2 >= 0) {
                outVec[0] = quartet[0]; outVec[1] = quartet[2];
                outVec[patternsI+1] = numToString(ABBAtotal); outVec[patternsI+2] = numToString(BBAAtotal);
                if (fStats) {
                    double Dnum = ABBAtotal - BBAAtotal;
                    outVec[7] = numToString(Dnum/F_G_denom2);
                }
            } else {
                outVec[0] = quartet[2]; outVec[1] = quartet[0];
                outVec[patternsI+1] = numToString(BBAAtotal); outVec[patternsI+2] = numToString(ABBAtotal);
                if (fStats) {
                    double Dnum = BBAAtotal - ABBAtotal;
                    outVec[7] = numToString(Dnum/F_G_denom2_reversed);
                }
            } break;
                
        case P3isTrios0:
            outVec[2] = quartet[0]; outVec[3] = quartet[3];
            outVec[4] = numToString(std::fabs(D3)); outVec[5] = numToString(D3_Z);
            outVec[6] = numToString(D3_p); outVec[patternsI] = numToString(ABBAtotal);
            if (D3 >= 0) {
                outVec[0] = quartet[2]; outVec[1] = quartet[1];
                outVec[patternsI+1] = numToString(BBAAtotal); outVec[patternsI+2] = numToString(BABAtotal);
                if (fStats) {
                    double Dnum = BBAAtotal - BABAtotal;
                    outVec[7] = numToString(Dnum/F_G_denom3);
                }
            } else {
                outVec[0] = quartet[1]; outVec[1] = quartet[2];
                outVec[patternsI+1] = numToString(BABAtotal); outVec[patternsI+2] = numToString(BBAAtotal);
                if (fStats) {
                    double Dnum = BABAtotal - BBAAtotal;
                    outVec[7] = numToString(Dnum/F_G_denom3_reversed);
                }
            } break;
                
        }
        return outVec;
    }
    
};


#endif /* Dsuite_utils_h */

