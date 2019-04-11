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
#include "gzstream.h"

#define PROGRAM_BIN "evo"
#define PACKAGE_BUGREPORT "mm21@sanger.ac.uk"
#define GZIP_EXT ".gz"

using std::string;
// VCF format constant
static const int NUM_NON_GENOTYPE_COLUMNS=9;  // 8 mendatory columns + 1 column with definition of the genotype columns


std::string stripExtension(const std::string& filename);
std::vector<std::string> split(const std::string &s, char delim);
std::vector<size_t> locateSet(std::vector<std::string>& sample_names, const std::vector<std::string>& set);
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode);
std::ostream* createWriter(const std::string& filename, std::ios_base::openmode mode);

// Converting numbers (int, double, size_t, and char) to string
template <typename T> std::string numToString(T i) {
    std::string ret;
    std::stringstream out;
    out << i;
    ret = out.str();
    return ret;
}

// Print an arbitrary vector to an output stream
template <class T> void print_vector_stream(T vector, std::ostream& outStream, char delim = '\t') {
    for (int i = 0; i < vector.size(); i++) {
        if (i == (vector.size()-1))
            outStream << vector[i] << std::endl;
        else
            outStream << vector[i] << delim;
    }
}
#endif /* Dsuite_utils_h */

