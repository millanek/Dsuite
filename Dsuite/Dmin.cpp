//
//  Dmin.cpp
//  Dsuite
//
//  Created by Milan Malinsky on 02/04/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "Dmin.h"

#define SUBPROGRAM "Dtrios"

#define DEBUG 0

static const char *DMIN_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf SETS.txt\n"
"Calculate the Dmin-statistic - the ABBA/BABA stat for all trios of species in the dataset (the outgroup being fixed)\n"
"the calculation is as definded in Durand et al. 2011\n"
"The SETS.txt should have two columns: SAMPLE_ID    SPECIES_ID\n"
"The outgroup (can be multiple samples) should be specified by using the keywork Outgroup in place of the SPECIES_ID\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -r , --region=start,length              (optional) only process a subset of the VCF file\n"
"       -t , --tree=TREE_FILE.nwk               (optional) a file with a tree in the newick format specifying the relationships between populations/species\n"
"                                               D values for trios arranged according to these relationships will be output in a file with _tree.txt suffix\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


static const char* shortopts = "hr:n:t:";

static const int JK_WINDOW = 20000;

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "region",   no_argument, NULL, 'r' },
    { "tree",   required_argument, NULL, 't' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string treeFile = "";
    static string runName = "";
    static int windowSize = 50;
    int jkWindowSize = JK_WINDOW;
    int regionStart = -1;
    int regionLength = -1;
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

int DminMain(int argc, char** argv) {
    parseDminOptions(argc, argv);
    string line; // for reading the input files
    std::istream* treeFile;
    std::ofstream* outFileTree;
    std::map<string,std::vector<int>> treeTaxonNamesToLoc; std::vector<int> treeLevels;
    if (opt::treeFile != "") {
        outFileTree = new std::ofstream(setsFileRoot+ "_" + opt::runName + "_tree.txt");
        treeFile = new std::ifstream(opt::treeFile.c_str());
        //std::regex branchLengths(":[0-9].[0-9]+");
        getline(*treeFile, line);
        treeLevels.assign(line.length(),0); int currentLevel = 0;
        std::vector<string> treeTaxonNames;
        string currentTaxonName = "";
        int lastBegin = 0;
        for (int i = 0; i < line.length(); ++i) {
            if (line[i] == '(') {
                currentLevel++; treeLevels[i] = currentLevel;
            } else if (line[i] == ')') {
                currentLevel--; treeLevels[i] = currentLevel;
                if (currentTaxonName != "") {
                    treeTaxonNames.push_back(currentTaxonName);
                    treeTaxonNamesToLoc[currentTaxonName].push_back(lastBegin);
                    treeTaxonNamesToLoc[currentTaxonName].push_back(i-1);
                    currentTaxonName = "";
                }
            } else if (line[i] == ',') {
                treeLevels[i] = currentLevel;
                if (currentTaxonName != "") {
                    treeTaxonNames.push_back(currentTaxonName);
                    treeTaxonNamesToLoc[currentTaxonName].push_back(lastBegin);
                    treeTaxonNamesToLoc[currentTaxonName].push_back(i-1);
                    currentTaxonName = "";
                }
            } else {
                if (currentTaxonName == "")
                    lastBegin = i;
                treeLevels[i] = currentLevel;
                currentTaxonName += line[i];
            }
        }
        //print_vector(treeTaxonNames, std::cout,'\n');
        //print_vector(treeLevels, std::cout,' ');
        //for (std::map<string,std::vector<int>>::iterator i = treeTaxonNamesToLoc.begin(); i != treeTaxonNamesToLoc.end(); i++) {
        //    std::cout << i->first << "\t" << i->second[0] << "\t" << i->second[1] << "\t" << treeLevels[i->second[0]] << "\t" << treeLevels[i->second[1]] << std::endl;
        //}
    }
    
    std::istream* vcfFile = createReader(opt::vcfFile.c_str());
    std::ifstream* setsFile = new std::ifstream(opt::setsFile.c_str());
    string setsFileRoot = stripExtension(opt::setsFile);
    std::ofstream* outFileBBAA;
    std::ofstream* outFileDmin;
    std::ofstream* outFileCombine;
    std::ofstream* outFileCombineStdErr;
    if (opt::regionStart == -1) {
        outFileBBAA = new std::ofstream(setsFileRoot+ "_" + opt::runName + "_BBAA.txt");
        outFileDmin = new std::ofstream(setsFileRoot+ "_" + opt::runName + "_Dmin.txt");
        outFileCombine = new std::ofstream(setsFileRoot+ "_" + opt::runName + "_combine.txt");
        outFileCombineStdErr = new std::ofstream(setsFileRoot+ "_" + opt::runName + "_combine_stderr.txt");
    } else {
        string fileNameString = setsFileRoot+"_"+opt::runName+"_"+numToString(opt::regionStart)+"_"+numToString(opt::regionStart+opt::regionLength);
        outFileBBAA = new std::ofstream(fileNameString+"_BBAA.txt");
        outFileDmin = new std::ofstream(fileNameString+"_Dmin.txt");
        outFileCombine = new std::ofstream(fileNameString+"_combine.txt");
        outFileCombineStdErr = new std::ofstream(fileNameString+"_combine_stderr.txt");
    }

    
    std::map<string, std::vector<string>> speciesToIDsMap;
    std::map<string, string> IDsToSpeciesMap;
    std::map<string, std::vector<size_t>> speciesToPosMap;
    std::map<size_t, string> posToSpeciesMap;
    
    // Get the sample sets
    bool outgroupSpecified = false;
    while (getline(*setsFile, line)) {
        // std::cerr << line << std::endl;
        std::vector<string> ID_Species = split(line, '\t');
        if (ID_Species[1] == "Outgroup") { outgroupSpecified = true; }
        speciesToIDsMap[ID_Species[1]].push_back(ID_Species[0]);
        IDsToSpeciesMap[ID_Species[0]] = ID_Species[1];
        //std::cerr << ID_Species[1] << "\t" << ID_Species[0] << std::endl;
    }
    if (!outgroupSpecified) { std::cerr << "The file " << opt::setsFile << " needs to specify the \"Outgroup\"" << std::endl; exit(1); }
    
    // Get a vector of set names (usually species)
    std::vector<string> species;
    for(std::map<string,std::vector<string>>::iterator it = speciesToIDsMap.begin(); it != speciesToIDsMap.end(); ++it) {
        if ((it->first) != "Outgroup" && it->first != "xxx") {
            species.push_back(it->first);
            // std::cerr << it->first << std::endl;
        }
    } std::cerr << "There are " << species.size() << " sets (excluding the Outgroup)" << std::endl;
    int nCombinations = nChoosek((int)species.size(),3);
    std::cerr << "Going to calculate " << nCombinations << " Dmin values" << std::endl;
    if (opt::treeFile != "") {
        
    }
    
    
    // first, get all combinations of three sets (species):
    std::vector<std::vector<string>> trios; trios.resize(nCombinations);
    std::vector<std::vector<int>> triosInt; triosInt.resize(nCombinations);
    std::vector<bool> v(species.size()); std::fill(v.begin(), v.begin() + 3, true); // prepare a selection vector
    int pNum = 0;
    do {
        for (int i = 0; i < v.size(); ++i) {
            if (v[i]) { trios[pNum].push_back(species[i]); triosInt[pNum].push_back(i); }
        } pNum++;
    } while (std::prev_permutation(v.begin(), v.end())); // Getting all permutations of the selection vector - so it selects all combinations
    std::cerr << "Done permutations" << std::endl;
    
    // And need to prepare the vectors to hold the D values:
    std::vector<double> init(3,0.0); // Vector of initial values
    std::vector<std::vector<double>> initDs(3); // vector with three empty (double) vectors
    std::vector<std::vector<std::vector<double>>> regionDs; regionDs.assign(nCombinations, initDs);
    std::vector<double> allPs(species.size(),0.0);
    std::vector<double> ABBAtotals(nCombinations,0); std::vector<double> BABAtotals(nCombinations,0);
    std::vector<double> BBAAtotals(nCombinations,0);
    std::vector<double> localABBAtotals(nCombinations,0); std::vector<double> localBABAtotals(nCombinations,0);
    std::vector<double> localBBAAtotals(nCombinations,0);
    std::vector<int> usedVars(nCombinations,0); // Will count the number of used variants for each trio
    int totalVariantNumber = 0;
    std::vector<string> sampleNames; std::vector<std::string> fields;
    // Find out how often to report progress, based on the number of trios
    int reportProgressEvery; if (nCombinations < 1000) reportProgressEvery = 100000;
    else if (nCombinations < 100000) reportProgressEvery = 10000;
    else reportProgressEvery = 1000;
    clock_t start; clock_t startGettingCounts; clock_t startCalculation;
    double durationOverall; double durationGettingCounts; double durationCalculation;
    
    while (getline(*vcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            // print_vector_stream(sampleNames, std::cerr);
            for (std::vector<std::string>::size_type i = 0; i != sampleNames.size(); i++) {
                posToSpeciesMap[i] = IDsToSpeciesMap[sampleNames[i]];
            }
            // Iterate over all the keys in the map to find the samples in the VCF:
            // Give an error if no sample is found for a species:
            for(std::map<string, std::vector<string>>::iterator it = speciesToIDsMap.begin(); it != speciesToIDsMap.end(); ++it) {
                string sp =  it->first;
                //std::cerr << "sp " << sp << std::endl;
                std::vector<string> IDs = it->second;
                std::vector<size_t> spPos = locateSet(sampleNames, IDs);
                if (spPos.empty()) {
                    std::cerr << "Did not find any samples in the VCF for \"" << sp << "\"" << std::endl;
                    assert(!spPos.empty());
                }
                speciesToPosMap[sp] = spPos;
            }
            start = clock();
            //  std::cerr << " " << std::endl;
            //  std::cerr << "Outgroup at pos: "; print_vector_stream(speciesToPosMap["Outgroup"], std::cerr);
        } else {
            totalVariantNumber++;
            if (opt::regionStart != -1) {
                if (totalVariantNumber < opt::regionStart)
                    continue;
                if (totalVariantNumber > (opt::regionStart+opt::regionLength)) {
                    std::cerr << "DONE" << std::endl; break;
                }
            }
            if (totalVariantNumber % reportProgressEvery == 0) {
                durationOverall = ( clock() - start ) / (double) CLOCKS_PER_SEC;
                std::cerr << "Processed " << totalVariantNumber << " variants in " << durationOverall << "secs" << std::endl;
                std::cerr << "GettingCounts " << durationGettingCounts << " calculation " << durationCalculation << "secs" << std::endl;
            }
            fields = split(line, '\t');
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());

            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4];
            if (refAllele.length() > 1 || altAllele.length() > 1) {
                refAllele.clear(); refAllele.shrink_to_fit(); altAllele.clear(); altAllele.shrink_to_fit();
                genotypes.clear(); genotypes.shrink_to_fit(); continue;
            }
            
            startGettingCounts = clock();
            GeneralSetCounts* c = new GeneralSetCounts(speciesToPosMap, (int)genotypes.size());
            c->getSetVariantCounts(genotypes, posToSpeciesMap);
            genotypes.clear(); genotypes.shrink_to_fit();
            durationGettingCounts = ( clock() - startGettingCounts ) / (double) CLOCKS_PER_SEC;
            
            startCalculation = clock();
            double p_O = c->setDAFs.at("Outgroup");
            if (p_O == -1) { delete c; continue; } // We need to make sure that the outgroup is defined
            
            for (std::vector<std::string>::size_type i = 0; i != species.size(); i++) {
                allPs[i] = c->setDAFs.at(species[i]);
            }
            
            // Now calculate the D stats:
            double p_S1; double p_S2; double p_S3; double ABBA; double BABA; double BBAA;
            for (int i = 0; i != trios.size(); i++) {
                p_S1 = allPs[triosInt[i][0]];
                if (p_S1 == -1) continue;  // If any member of the trio has entirely missing data, just move on to the next trio
                p_S2 = allPs[triosInt[i][1]];
                if (p_S2 == -1) continue;
                p_S3 = allPs[triosInt[i][2]];
                if (p_S3 == -1) continue;
                usedVars[i]++;
                
                ABBA = ((1-p_S1)*p_S2*p_S3*(1-p_O)); ABBAtotals[i] += ABBA; localABBAtotals[i] += ABBA;
                BABA = (p_S1*(1-p_S2)*p_S3*(1-p_O)); BABAtotals[i] += BABA; localBABAtotals[i] += BABA;
                BBAA = ((1-p_S3)*p_S2*p_S1*(1-p_O)); BBAAtotals[i] += BBAA; localBBAAtotals[i] += BBAA;
                
                

                
                if (usedVars[i] % opt::jkWindowSize == 0) {
                    double localDnums1 = localABBAtotals[i] - localBABAtotals[i]; double localDnums2 = localABBAtotals[i] - localBBAAtotals[i]; double localDnums3 = localBBAAtotals[i] - localBABAtotals[i];
                    double localDdenoms1 = localABBAtotals[i] + localBABAtotals[i]; double localDdenoms2 = localABBAtotals[i] + localBBAAtotals[i]; double localDdenoms3 = localBBAAtotals[i] + localBABAtotals[i];
                    double regionD0 = localDnums1/localDdenoms1; double regionD1 = localDnums2/localDdenoms2;
                    double regionD2 = localDnums3/localDdenoms3;
                    regionDs[i][0].push_back(regionD0); regionDs[i][1].push_back(regionD1); regionDs[i][2].push_back(regionD2);
                    localABBAtotals[i] = 0; localBABAtotals[i] = 0; localBBAAtotals[i] = 0;
                }
                // }
            }
            durationCalculation = ( clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
            delete c;
        }
    }
    
    
    for (int i = 0; i != trios.size(); i++) { //
        // Get the standard error values:
        double D1stdErr = jackknive_std_err(regionDs[i][0]); double D2stdErr = jackknive_std_err(regionDs[i][1]);
        double D3stdErr = jackknive_std_err(regionDs[i][2]);
        // Get the D values
        double Dnum1 = ABBAtotals[i] - BABAtotals[i];
        double Dnum2 = ABBAtotals[i] - BBAAtotals[i];
        double Dnum3 = BBAAtotals[i] - BABAtotals[i];
        
        double Ddenom1 = ABBAtotals[i] + BABAtotals[i];
        double Ddenom2 = ABBAtotals[i] + BBAAtotals[i];
        double Ddenom3 = BBAAtotals[i] + BABAtotals[i];
        double D1 = Dnum1/Ddenom1; double D2 = Dnum2/Ddenom2; double D3 = Dnum3/Ddenom3;
        // Get the Z-scores
        double D1_Z = fabs(D1)/D1stdErr; double D2_Z = fabs(D2)/D2stdErr;
        double D3_Z = fabs(D3)/D3stdErr;
        
        
        // Find which topology is in agreement with the counts of the BBAA, BABA, and ABBA patterns
        if (BBAAtotals[i] >= BABAtotals[i] && BBAAtotals[i] >= ABBAtotals[i]) {
            if (D1 >= 0)
                *outFileBBAA << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2];
            else
                *outFileBBAA << trios[i][1] << "\t" << trios[i][0] << "\t" << trios[i][2];
            *outFileBBAA << "\t" << fabs(D1) << "\t" << D1_Z << "\t";
            *outFileBBAA << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
        } else if (BABAtotals[i] >= BBAAtotals[i] && BABAtotals[i] >= ABBAtotals[i]) {
            if (D2 >= 0)
                *outFileBBAA << trios[i][0] << "\t" << trios[i][2] << "\t" << trios[i][1];
            else
                *outFileBBAA << trios[i][2] << "\t" << trios[i][0] << "\t" << trios[i][1];
            *outFileBBAA << "\t" << fabs(D2) << "\t" << D2_Z << "\t";
            *outFileBBAA << BABAtotals[i] << "\t" << BBAAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
        } else if (ABBAtotals[i] >= BBAAtotals[i] && ABBAtotals[i] >= BABAtotals[i]) {
            if (D3 >= 0)
                *outFileBBAA << trios[i][2] << "\t" << trios[i][1] << "\t" << trios[i][0];
            else
                *outFileBBAA << trios[i][1] << "\t" << trios[i][2] << "\t" << trios[i][0];
            *outFileBBAA << "\t" << fabs(D3) << "\t" << D3_Z << "\t";
            *outFileBBAA << ABBAtotals[i] << "\t" << BABAtotals[i] << "\t" << BBAAtotals[i] << std::endl;
        }
        
        // Find Dmin:
        if (fabs(D1) <= fabs(D2) && fabs(D1) <= fabs(D3)) { // (P3 == S3)
            if (D1 >= 0)
                *outFileDmin << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << D1 << "\t" << D1_Z << "\t" << std::endl;
            else
                *outFileDmin << trios[i][1] << "\t" << trios[i][0] << "\t" << trios[i][2] << "\t" << fabs(D1) << "\t" << D1_Z << "\t"<< std::endl;
            // if (BBAAtotals[i] < BABAtotals[i] || BBAAtotals[i] < ABBAtotals[i])
            //     std::cerr << "\t" << "WARNING: Dmin tree different from DAF tree" << std::endl;
        } else if (fabs(D2) <= fabs(D1) && fabs(D2) <= fabs(D3)) { // (P3 == S2)
            if (D2 >= 0)
                *outFileDmin << trios[i][0] << "\t" << trios[i][2] << "\t" << trios[i][1] << "\t" << D2 << "\t" << D2_Z << "\t"<< std::endl;
            else
                *outFileDmin << trios[i][2] << "\t" << trios[i][0] << "\t" << trios[i][1] << "\t" << fabs(D2) << "\t" << D2_Z << "\t"<< std::endl;
            // if (BABAtotals[i] < BBAAtotals[i] || BABAtotals[i] < ABBAtotals[i])
            //     std::cerr << "\t" << "WARNING: Dmin tree different from DAF tree" << std::endl;
        } else if (fabs(D3) <= fabs(D1) && fabs(D3) <= fabs(D2)) { // (P3 == S1)
            if (D3 >= 0)
                *outFileDmin << trios[i][2] << "\t" << trios[i][1] << "\t" << trios[i][0] << "\t" << D3 << "\t" << D3_Z << "\t"<< std::endl;
            else
                *outFileDmin << trios[i][1] << "\t" << trios[i][2] << "\t" << trios[i][0] << "\t" << fabs(D3) << "\t" << D3_Z << "\t" << std::endl;
            // if (ABBAtotals[i] < BBAAtotals[i] || ABBAtotals[i] < BABAtotals[i])
            //     std::cerr << "\t" << "WARNING: Dmin tree different from DAF tree" << std::endl;
        }
        
        // Find which arrangement of trios is consistent with the input tree (if provided):
        if (opt::treeFile != "") {
            int loc1 = treeTaxonNamesToLoc[trios[i][0]][0];
            int loc2 = treeTaxonNamesToLoc[trios[i][1]][0];
            int loc3 = treeTaxonNamesToLoc[trios[i][2]][0];
            
            int arrangement = 0;    // 1 - trios[i][0] and trios[i][2] are P1 and P2
                                    // 2 - trios[i][0] and trios[i][1] are P1 and P2
                                    // 3 - trios[i][1] and trios[i][2] are P1 and P2
            int midLoc = std::max(std::min(loc1,loc2), std::min(std::max(loc1,loc2),loc3));
            if (midLoc == loc1) {
                if (loc2 < loc1) {
                    int m1 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc1);
                    int m2 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc3);
                    if (m1 < m2) arrangement = 1; else arrangement = 2;
                } else {
                    int m1 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc1);
                    int m2 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc2);
                    if (m1 < m2) arrangement = 2; else arrangement = 1;
                }
            } else if (midLoc == loc2) {
                if (loc1 < loc2) {
                    int m1 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc2);
                    int m2 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc3);
                    if (m1 < m2) arrangement = 3; else arrangement = 2;
                } else {
                    int m1 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc2);
                    int m2 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc1);
                    if (m1 < m2) arrangement = 2; else arrangement = 3;
                }
            } else if (midLoc == loc3) {
                if (loc1 < loc3) {
                    int m1 = *std::min_element(treeLevels.begin()+loc1, treeLevels.begin()+loc3);
                    int m2 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc2);
                    if (m1 < m2) arrangement = 3; else arrangement = 1;
                    
                } else {
                    int m1 = *std::min_element(treeLevels.begin()+loc2, treeLevels.begin()+loc3);
                    int m2 = *std::min_element(treeLevels.begin()+loc3, treeLevels.begin()+loc1);
                    if (m1 < m2) arrangement = 1; else arrangement = 3;
                }
            }
            switch (arrangement)
            {
                case 1:
                    if (D2 >= 0)
                        *outFileTree << trios[i][0] << "\t" << trios[i][2] << "\t" << trios[i][1] << "\t" << D2 << "\t" << D2_Z << "\t"<< std::endl;
                    else
                        *outFileTree << trios[i][2] << "\t" << trios[i][0] << "\t" << trios[i][1] << "\t" << fabs(D2) << "\t" << D2_Z << "\t"<< std::endl;
                    break;
                case 2:
                    if (D1 >= 0)
                        *outFileTree << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << D1 << "\t" << D1_Z << "\t" << std::endl;
                    else
                        *outFileTree << trios[i][1] << "\t" << trios[i][0] << "\t" << trios[i][2] << "\t" << fabs(D1) << "\t" << D1_Z << "\t"<< std::endl;
                    break;
                case 3:
                    if (D3 >= 0)
                        *outFileTree << trios[i][2] << "\t" << trios[i][1] << "\t" << trios[i][0] << "\t" << D3 << "\t" << D3_Z << "\t"<< std::endl;
                    else
                        *outFileTree << trios[i][1] << "\t" << trios[i][2] << "\t" << trios[i][0] << "\t" << fabs(D3) << "\t" << D3_Z << "\t" << std::endl;
                    break;
            }

        }
        
        // Output a simple file that can be used for combining multiple local runs:
        *outFileCombine << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
        print_vector(regionDs[i][0], *outFileCombineStdErr, ',', false); *outFileCombineStdErr << "\t"; print_vector(regionDs[i][1], *outFileCombineStdErr, ',', false); *outFileCombineStdErr << "\t";
        print_vector(regionDs[i][2], *outFileCombineStdErr, ',',false); *outFileCombineStdErr << std::endl;
        
        //std::cerr << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << D1 << "\t" << D2 << "\t" << D3 << "\t" << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
    }
    return 0;
    
}



void parseDminOptions(int argc, char** argv) {
    bool die = false; string regionArgString; std::vector<string> regionArgs;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 't': arg >> opt::treeFile; break;
            case 'r': arg >> regionArgString; regionArgs = split(regionArgString, ',');
                opt::regionStart = (int)stringToDouble(regionArgs[0]); opt::regionLength = (int)stringToDouble(regionArgs[1]);  break;
            case 'h':
                std::cout << DMIN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 2) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 2)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DMIN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
}

