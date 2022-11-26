//
//  D.cpp
//  Dsuite
//
//  Created by Milan Malinsky on 11/04/2019.
//

#include "D.h"
#include "Dsuite_common.h"
#include "kstest.h"
#include <deque>
#include <list>
#define SUBPROGRAM "Dinvestigate"

#define DEBUG 1
#define MIN_SETS 3

static const char *ABBA_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf.gz SETS.txt test_trios.txt\n"
"Outputs D, f_d (Martin et al. 2014 MBE), f_dM (Malinsky et al., 2015), and d_f (Pfeifer & Kapan, 2019) in genomic windows\n"
"The SETS.txt file should have two columns: SAMPLE_ID    POPULATION_ID\n"
"The test_trios.txt should contain names of three populations for which the statistics will be calculated:\n"
"POP1   POP2    POP3\n"
"There can be multiple lines and then the program generates multiple ouput files, named like POP1_POP2_POP3_localFstats_SIZE_STEP.txt\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -w SIZE,STEP --window=SIZE,STEP         (required) D, f_D, f_dM, and d_f statistics for windows containing SIZE useable SNPs, moving by STEP (default: 50,25)\n"
"       -g, --use-genotype-probabilities        (optional) use probabilities (GP tag) or calculate them from likelihoods (GL or PL tags) using a Hardy-Weinberg prior\n"
"                                               the probabilities are used to estimate allele frequencies in each population/species\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


//enum { OPT_F_JK };

static const char* shortopts = "hw:n:g";

//static const int JK_WINDOW = 5000;

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "window",   required_argument, NULL, 'w' },
    { "help",   no_argument, NULL, 'h' },
    { "use-genotype-probabilities", no_argument, NULL, 'g'},
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string testTriosFile;
    static string runName = "";
    static int minScLength = 0;
    static int windowSize = 50;
    static int windowStep = 25;
    static bool useGenotypeProbabilities = false;
    //int jkWindowSize = JK_WINDOW;
}


void doAbbaBaba() {
    string line; // for reading the input files
    
    std::istream* vcfFile = createReader(opt::vcfFile);
    std::ifstream* testTriosFile = new std::ifstream(opt::testTriosFile.c_str());
    if (!testTriosFile->good()) { std::cerr << "The file " << opt::testTriosFile << " could not be opened. Exiting..." << std::endl; exit(EXIT_FAILURE);}
    
    // Get the sample sets
    SetInformation setInfo(opt::setsFile, MIN_SETS, OutgroupRequired);
    
    // Get the test trios
    std::vector<std::ofstream*> outFiles;
    std::vector<std::ofstream*> outFilesGenes;
    std::vector<std::vector<string> > testTrios;
    while (getline(*testTriosFile,line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        // std::cerr << line << std::endl;
        std::vector<string> threePops = split(line, '\t'); assert(threePops.size() == 3);
        for (int i = 0; i != threePops.size(); i++) { // Check that the test trios are in the sets file
            if (setInfo.popToIDsMap.count(threePops[i]) == 0) {
                std::cerr << threePops[i] << " is present in the " << opt::testTriosFile << " but missing from the " << opt::setsFile << std::endl;
            }
        }
        std::ofstream* outFile = new std::ofstream(threePops[0] + "_" + threePops[1] + "_" + threePops[2]+ "_localFstats_" + opt::runName + "_" + numToString(opt::windowSize) + "_" + numToString(opt::windowStep) + ".txt");
        *outFile << "chr\twindowStart\twindowEnd\tD\tf_d\tf_dM\td_f" << std::endl;
        outFiles.push_back(outFile);
        testTrios.push_back(threePops);
    }
    
    // Create objects to hold the results for each trio
    TestTrioInfo info(opt::windowSize); std::vector<TestTrioInfo> testTrioInfos(testTrios.size(), info);
    
    // Now go through the vcf and calculate D
    int totalVariantNumber = 0;
    int reportProgressEvery = 1000; string chr; string coord;

   // int lastPrint = 0; int lastWindowVariant = 0;
    std::vector<string> sampleNames; std::vector<std::string> fields;
    clock_t start; clock_t startGettingCounts; clock_t startCalculation;
    double durationOverall; double durationGettingCounts; double durationCalculation;
    while (getline(*vcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            setInfo.linkSetsAndVCFpositions(sampleNames);
            start = clock();
        } else {
            totalVariantNumber++;
            if (totalVariantNumber % reportProgressEvery == 0) {
                durationOverall = ( clock() - start ) / (double) CLOCKS_PER_SEC;
                std::cerr << "Processed " << totalVariantNumber << " variants in " << durationOverall << "secs" << std::endl;
                //std::cerr << "GettingCounts " << durationGettingCounts << " calculation " << durationCalculation << "secs" << std::endl;
            }
            fields = split(line, '\t'); chr = fields[0]; coord = fields[1];
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4];
            if (refAllele.length() > 1 || altAllele.length() > 1 || altAllele == "*") {
                refAllele.clear(); refAllele.shrink_to_fit(); altAllele.clear(); altAllele.shrink_to_fit();
                genotypes.clear(); genotypes.shrink_to_fit(); continue;
            }
            
            startGettingCounts = clock();
            GeneralSetCounts* c = new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size());
            try { c->getSetVariantCounts(genotypes, setInfo.posToPopMap); } catch (const std::out_of_range& oor) {
                std::cerr << "Problems getting splitCounts for " << chr << " " << coord << std::endl; }
            if (opt::useGenotypeProbabilities) {
                int likelihoodsOrProbabilitiesTagPosition = c->checkForGenotypeLikelihoodsOrProbabilities(fields);
                if (likelihoodsOrProbabilitiesTagPosition == LikelihoodsProbabilitiesAbsent) {
                    printMissingLikelihoodsWarning(fields[0], fields[1]);
                    opt::useGenotypeProbabilities = false;
                } else c->getAFsFromGenotypeLikelihoodsOrProbabilities(genotypes,setInfo.posToPopMap,likelihoodsOrProbabilitiesTagPosition);
            }
            genotypes.clear(); genotypes.shrink_to_fit();
            durationGettingCounts = ( clock() - startGettingCounts ) / (double) CLOCKS_PER_SEC;
            
            startCalculation = clock();
            double p_O; try { p_O = c->setDAFs.at("Outgroup"); } catch (const std::out_of_range& oor) {
                std::cerr << "Counts don't contain derived allele frequency for the Outgroup" << std::endl; }
            if (p_O == -1) { delete c; continue; } // We need to make sure that the outgroup is defined
            
            double p_S1; double p_S2; double p_S3; double ABBA; double BABA; double F_d_denom; double F_dM_denom;
            for (int i = 0; i != testTrios.size(); i++) {
                try {
                    if (!opt::useGenotypeProbabilities) p_S1 = c->setDAFs.at(testTrios[i][0]);
                    else p_S1 = c->setDAFsFromLikelihoods.at(testTrios[i][0]);
                } catch (const std::out_of_range& oor) {
                std::cerr << "Counts don't contain derived allele frequency for " << testTrios[i][0] << std::endl; }
                if (p_S1 == -1) continue;  // If any member of the trio has entirely missing data, just move on to the next trio
                try {
                    if (!opt::useGenotypeProbabilities) p_S2 = c->setDAFs.at(testTrios[i][1]);
                    else p_S2 = c->setDAFsFromLikelihoods.at(testTrios[i][1]);
                } catch (const std::out_of_range& oor) {
                    std::cerr << "Counts don't contain derived allele frequency for " << testTrios[i][1] << std::endl; }
                if (p_S2 == -1) continue;
                try {
                    if (!opt::useGenotypeProbabilities) p_S3 = c->setDAFs.at(testTrios[i][2]);
                    else p_S3 = c->setDAFsFromLikelihoods.at(testTrios[i][2]);
                } catch (const std::out_of_range& oor) {
                    std::cerr << "Counts don't contain derived allele frequency for " << testTrios[i][2] << std::endl; }
                if (p_S3 == -1) continue;
                //if (p_S3 == 0) continue; // XXAA pattern is not informative
                if (p_S1 == 0 && p_S2 == 0 && p_S3 == 0) continue; // Checking if the SNP is variable in the trio
                if (p_S1 == 1 && p_S2 == 1 && p_S3 == 1) continue; // Checking if the SNP is variable in the trio
                //if (p_S1 == 1 && p_S2 == 1) continue; // BBAA pattern is not informative
                //if (p_S1 == 0 && p_S2 == 0) continue; // AABA pattern is not informative
                
                
                ABBA = ((1-p_S1)*p_S2*p_S3*(1-p_O)); testTrioInfos[i].ABBAtotal += ABBA;
                if(ABBA > 0.5) {
                    testTrioInfos[i].ABBAsitePositionsPerChomosome[chr].push_back(atoi(coord.c_str()));
                }
                BABA = (p_S1*(1-p_S2)*p_S3*(1-p_O)); testTrioInfos[i].BABAtotal += BABA;
                if(BABA > 0.5) {
                    testTrioInfos[i].BABAsitePositionsPerChomosome[chr].push_back(atoi(coord.c_str()));
                }
                
                if (p_S2 > p_S3) {
                    F_d_denom = ((1-p_S1)*p_S2*p_S2*(1-p_O)) - (p_S1*(1-p_S2)*p_S2*(1-p_O));
                } else {
                    F_d_denom = ((1-p_S1)*p_S3*p_S3*(1-p_O)) - (p_S1*(1-p_S3)*p_S3*(1-p_O));
                } testTrioInfos[i].F_d_denom += F_d_denom; testTrioInfos[i].interimF_d_denom += F_d_denom;
                
                if (p_S1 <= p_S2) {
                    if (p_S2 > p_S3) {
                        F_dM_denom = ((1-p_S1)*p_S2*p_S2*(1-p_O)) - (p_S1*(1-p_S2)*p_S2*(1-p_O));
                    } else {
                        F_dM_denom = ((1-p_S1)*p_S3*p_S3*(1-p_O)) - (p_S1*(1-p_S3)*p_S3*(1-p_O));
                    }
                } else {
                    if (p_S1 > p_S3) {
                        F_dM_denom = -(((1-p_S1)*p_S2*p_S1*(1-p_O)) - (p_S1*(1-p_S2)*p_S1)*(1-p_O));
                    } else {
                        F_dM_denom = -(((1-p_S3)*p_S2*p_S3*(1-p_O)) - (p_S3*(1-p_S2)*p_S3)*(1-p_O));
                    }
                } testTrioInfos[i].F_dM_denom += F_dM_denom; testTrioInfos[i].interimF_dM_denom += F_dM_denom;
                
                
                // d_f
                double d13 = p_S1 + p_S3 - (2*p_S1*p_S3); double d23 = p_S2 + p_S3 - (2*p_S2*p_S3);
                double dfNum = p_S2 * d13 - p_S1 * d23;
                double dfDenom = p_S2 * d13 + p_S1 * d23;
                
                double ABBAplusBABA = ABBA + BABA;
                if (ABBAplusBABA != 0) {
                    testTrioInfos[i].windowABBAs.push_back(ABBA);  testTrioInfos[i].windowBABAs.push_back(BABA);
                    testTrioInfos[i].windowF_d_denoms.push_back(testTrioInfos[i].interimF_d_denom);
                    testTrioInfos[i].windowF_dM_denoms.push_back(testTrioInfos[i].interimF_dM_denom);
                    testTrioInfos[i].window_d_f_nums.push_back(dfNum); testTrioInfos[i].window_d_f_denoms.push_back(dfDenom);
                    testTrioInfos[i].windowInformativeSitesCords.push_back(atoi(coord.c_str()));
                    testTrioInfos[i].windowABBAs.pop_front(); testTrioInfos[i].windowBABAs.pop_front();
                    testTrioInfos[i].windowF_d_denoms.pop_front(); testTrioInfos[i].windowF_dM_denoms.pop_front();
                    testTrioInfos[i].windowInformativeSitesCords.pop_front();
                    testTrioInfos[i].window_d_f_nums.pop_front(); testTrioInfos[i].window_d_f_denoms.pop_front();
                    testTrioInfos[i].interimF_d_denom = 0; testTrioInfos[i].interimF_dM_denom = 0;
                    testTrioInfos[i].usedVars++;
                
                    if ((testTrioInfos[i].usedVars > opt::windowSize) && (testTrioInfos[i].usedVars % opt::windowStep == 0)) {
                        double windowABBAtotal = vector_sum(testTrioInfos[i].windowABBAs); double windowBABAtotal = vector_sum(testTrioInfos[i].windowBABAs);
                        double windowF_d_denom = vector_sum(testTrioInfos[i].windowF_d_denoms); double windowF_dM_denom = vector_sum(testTrioInfos[i].windowF_dM_denoms);
                        double wDnum = windowABBAtotal - windowBABAtotal; double wDdenom = windowABBAtotal + windowBABAtotal;
                        double w_d_f_num = vector_sum(testTrioInfos[i].window_d_f_nums);
                        double w_d_f_denom = vector_sum(testTrioInfos[i].window_d_f_denoms);
                        if ((atoi(coord.c_str()) - testTrioInfos[i].windowInformativeSitesCords[0]) > 0) {
                            *outFiles[i] << std::fixed << chr << "\t" << testTrioInfos[i].windowInformativeSitesCords[0] << "\t" << coord << "\t" << wDnum/wDdenom << "\t" << wDnum/windowF_d_denom << "\t" << wDnum/windowF_dM_denom << "\t" << w_d_f_num/w_d_f_denom << std::endl;
                        }
                    }
                }
            }
            durationCalculation = ( clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
            delete c;
        }
    }
    
    for (int i = 0; i != testTrios.size(); i++) {
        testTrioInfos[i].mergeABBA_BABA_SiteCoordsOverChoms(); testTrioInfos[i].testIfSitesUniformlyDistributed();
        
        std::cout << testTrios[i][0] << "\t" << testTrios[i][1] << "\t" << testTrios[i][2] << std::endl;
        std::cout << "D=" << (double)(testTrioInfos[i].ABBAtotal-testTrioInfos[i].BABAtotal)/(testTrioInfos[i].ABBAtotal+testTrioInfos[i].BABAtotal) << std::endl;
        std::cout << "f_d=" << (double)(testTrioInfos[i].ABBAtotal-testTrioInfos[i].BABAtotal)/testTrioInfos[i].F_d_denom << "\t" << (testTrioInfos[i].ABBAtotal-testTrioInfos[i].BABAtotal) << "/" << testTrioInfos[i].F_d_denom << std::endl;
        std::cout << "f_dM=" << (double)(testTrioInfos[i].ABBAtotal-testTrioInfos[i].BABAtotal)/testTrioInfos[i].F_dM_denom << "\t" << (testTrioInfos[i].ABBAtotal-testTrioInfos[i].BABAtotal) << "/" << testTrioInfos[i].F_dM_denom << std::endl;
        std::cout << "ABBA_KSpval = " << testTrioInfos[i].ABBA_KSpval << std::endl;
        std::cout << "BABA_KSpval = " << testTrioInfos[i].BABA_KSpval << std::endl;
        std::cout << std::endl;
    }
}


int abbaBabaMain(int argc, char** argv) {
    parseAbbaBabaOptions(argc, argv);
    doAbbaBaba();
    return 0;
    
}

void TestTrioInfo::testIfSitesUniformlyDistributed() {
    // Take care of the splits by random sampling with replacement:
    std::random_device rd;     // only used once to initialise (seed) engine
    std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
    std::uniform_int_distribution<int> uniABBA(0,linearABBApos.back()); // guaranteed unbiased
    std::uniform_int_distribution<int> uniBABA(0,linearBABApos.back()); // guaranteed unbiased
    std::list<int64_t> uniABBAvals; std::list<int64_t> uniBABAvals;
    // uniABBAvals.re(linearABBApos.size()); uniBABAvals.resize(linearBABApos.size());
    
    
    int numUniformSamples = (int)linearABBApos.size(); if (numUniformSamples < 10000) { numUniformSamples = 10000; }
    for (int i = 0; i < numUniformSamples; i++) {
        uniABBAvals.push_back(uniABBA(rng));
    }
    
    numUniformSamples = (int)linearBABApos.size(); if (numUniformSamples < 10000) { numUniformSamples = 10000; }
    for (int i = 0; i < numUniformSamples; i++) {
        uniBABAvals.push_back(uniBABA(rng));
    }
    
    std::list<int64_t> linearABBAposList(linearABBApos.begin(),linearABBApos.end());
    std::list<int64_t> linearBABAposList(linearBABApos.begin(),linearBABApos.end());
    
    ABBA_KSpval = ks_test(uniABBAvals, linearABBAposList, std::cerr, false);
    BABA_KSpval = ks_test(uniBABAvals, linearBABAposList, std::cerr, false);
    
    //double BABApval = ks_test(uniBABAvals, linearBABApos, std::cerr);
    
}

 


void TestTrioInfo::mergeABBA_BABA_SiteCoordsOverChoms() {
    int totalNumABBAsites = 0;
    for(std::map<string,std::vector<int>>::iterator it = ABBAsitePositionsPerChomosome.begin(); it != ABBAsitePositionsPerChomosome.end(); it++) {
        totalNumABBAsites = totalNumABBAsites + (int)it->second.size();
    } linearABBApos.reserve(totalNumABBAsites);
    
    int linearPosSoFar = 0;
    for(std::map<string,std::vector<int>>::iterator it = ABBAsitePositionsPerChomosome.begin(); it != ABBAsitePositionsPerChomosome.end(); it++) {
        for (std::vector<int>::size_type i = 0; i < it->second.size(); i++) {
            linearABBApos.push_back(it->second[i] + linearPosSoFar);
        }
    }
    
    int totalNumBABAsites = 0;
    for(std::map<string,std::vector<int>>::iterator it = BABAsitePositionsPerChomosome.begin(); it != BABAsitePositionsPerChomosome.end(); it++) {
        totalNumBABAsites = totalNumBABAsites + (int)it->second.size();
    } linearBABApos.reserve(totalNumBABAsites);
    
    linearPosSoFar = 0;
    for(std::map<string,std::vector<int>>::iterator it = BABAsitePositionsPerChomosome.begin(); it != BABAsitePositionsPerChomosome.end(); it++) {
        for (std::vector<int>::size_type i = 0; i < it->second.size(); i++) {
            linearBABApos.push_back(it->second[i] + linearPosSoFar);
        }
    }
    
}

void parseAbbaBabaOptions(int argc, char** argv) {
    bool die = false;
    std::vector<string> windowSizeStep;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'w':
                windowSizeStep = split(arg.str(), ',');
                if(windowSizeStep.size() != 2) {std::cerr << "The -w option requires two arguments, separated by a comma ','\n"; exit(EXIT_FAILURE);}
                opt::windowSize = atoi(windowSizeStep[0].c_str());
                opt::windowStep = atoi(windowSizeStep[1].c_str());
                break;
            case 'n': arg >> opt::runName; break;
            case 'g': opt::useGenotypeProbabilities = true; break;
            case 'h':
                std::cout << ABBA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (argc - optind < 3) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > 3)
    {
        std::cerr << "too many arguments\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << ABBA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
    opt::testTriosFile = argv[optind++];
}
