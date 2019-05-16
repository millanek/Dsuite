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
#include <math.h>
#include <algorithm>
#include <assert.h>
#include <time.h>
#include <regex>
#include "gzstream.h"

#define PROGRAM_BIN "Dsuite"
#define PACKAGE_BUGREPORT "milan.malinsky@unibas.ch"
#define GZIP_EXT ".gz"

using std::string;
// VCF format constant
static const int NUM_NON_GENOTYPE_COLUMNS=9;  // 8 mendatory columns + 1 column with definition of the genotype columns

double normalCDF(double x);
double stringToDouble(std::string s);
std::string stripExtension(const std::string& filename);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<size_t> locateSet(std::vector<std::string>& sample_names, const std::vector<std::string>& set);
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode = std::ios_base::in);
std::ostream* createWriter(const std::string& filename, std::ios_base::openmode mode = std::ios_base::out);
bool file_exists(const std::string& name);

// Converting numbers (int, double, size_t, and char) to string
template <typename T> std::string numToString(T i) {
    std::string ret;
    std::stringstream out;
    out << i;
    ret = out.str();
    return ret;
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
    if (vector.size() <= 2) {
        throw "Not enough blocks to calculate jackknife!!\n"
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
        sum += pow((jackkniveAverages[i] - jackkniveOverallMean), 2.0);
    }
    double var = ((double)(jackkniveAverages.size()-1)/(double)jackkniveAverages.size()) * sum;
    double Dstd_err = sqrt(var);
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
    void getBasicCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
};

#endif /* Dsuite_utils_h */

