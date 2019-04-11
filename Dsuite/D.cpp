//
//  D.cpp
//  Dsuite
//
//  Created by Milan Malinsky on 11/04/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#include "D.h"
#define SUBPROGRAM "abba-baba"

#define DEBUG 1

static const char *ABBA_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf SETS.txt\n"
"Calculate the D-statistic (abba/baba) as definded in Durand et al. 2011"
"This is a Four-Taxon Statistic to Test for Admixture using data from a VCF file\n"
"Now also calculates f_d as defined by Martin et al. 2014 MBE paper - for the ABBA pattern"
"The SETS.txt should have exactly four lines:\n"
"Line 1: Outgroup individuals\n"
"Line 2,3,4: P3,P2,P1 individuals respectively (as defined in the Durand et al. 2011 paper\n"
"\n"
"       -h, --help                              display this help and exit\n"
"       -f, --frequency                         use allele frequency data instead of single sequences for each of (P1,P2,P3,O)\n"
"       --AAeqO                                 ancestral allele info in the VCF is from the outgroup (e.g. Pnyererei for Malawi)\n"
"       --NoAaO                                 there is no ancestral allele info in the VCF AA field\n"
"       -w SIZE, --window=SIZE                  (optional) output D statistics for nonoverlapping windows containing SIZE SNPs with nonzero D (default: 50)\n"
"       -s SAMPLES.txt, --samples=SAMPLES.txt   (optional) supply a file of sample identifiers\n"
"                                               (default: sample ids from the vcf file are used)\n"
"       -n, --run-name                          run-name will be included in the output file name\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


enum { OPT_AA_EQ_O, OPT_NO_AA_O };

static const char* shortopts = "hs:fw:n:";

static const int JACKKNIVE_WINDOW_SIZE_FREQUENCY = 5000;
static const int JACKKNIVE_WINDOW_SIZE_SEQUENCE = 2000;

static const struct option longopts[] = {
    { "samples",   required_argument, NULL, 's' },
    { "run-name",   required_argument, NULL, 'n' },
    { "window",   required_argument, NULL, 'w' },
    { "AAeqO",   no_argument, NULL, OPT_AA_EQ_O },
    { "NoAaO",   no_argument, NULL, OPT_NO_AA_O },
    { "frequency",   no_argument, NULL, 'f' },
    { "help",   no_argument, NULL, 'h' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string sampleNameFile;
    static string runName = "";
    static bool bFrequency = false;
    static bool bAaEqO = false;
    static bool bNoAaO = false;
    static int minScLength = 0;
    static int windowSize = 50;
    int jackKniveWindowSize = JACKKNIVE_WINDOW_SIZE_SEQUENCE;
}

namespace ABBABABAcounts {
    int AABA = 0;
    int XXAA = 0;
    int BBBA = 0;
    int XXBA = 0;
    int ABBA = 0;
    int BABA = 0;
    int p1p2 = 0;
    int indels = 0;
    int noDafInfo = 0;
    int usedVariantsCounter = 0;
    int used_f_d_Counter = 0;
    // int ABBABABA = 0;
}


class ABBA_BABA_Freq_allResults {
public:
    ABBA_BABA_Freq_allResults() : Dnumerator(0), Ddenominator(0), lastVarsDnum(0), lastVarsDdenom(0), windowDnum(0), windowDdenom(0), f_d_denominator(0), window_f_d_denominator(0),  lastVarsF_d_denom(0), f_d_num(0), window_f_d_num(0),  lastVarsF_d_num(0),f_G_denom(0), f_G_num(0), lastVarsF_G_num(0), lastVarsF_G_denom(0), f_dM_denominator(0), window_f_dM_denominator(0), lastVarsF_dM_denom(0) {};
    double Dnumerator; double Ddenominator;     // simple D statistic
    double lastVarsDnum; double lastVarsDdenom; // D within a long stretch window for jackkinive analysis
    double windowDnum; double windowDdenom; // D within a window
    double f_d_denominator; double window_f_d_denominator; double lastVarsF_d_denom;
    double f_dM_denominator; double window_f_dM_denominator; double lastVarsF_dM_denom;
    double f_d_num; double window_f_d_num; double lastVarsF_d_num;
    double f_G_denom; double f_G_num; double lastVarsF_G_num; double lastVarsF_G_denom;
};

inline void incrementDnumDdenomFrequency(const ThreeSetCounts& c, ABBA_BABA_Freq_allResults& res) {
    if (c.set1daAF == -1) {
        ABBABABAcounts::noDafInfo++;
    } else if (c.set3daAF == 0) {
        ABBABABAcounts::XXAA++;
    } else if (c.set1daAF == 0 && c.set2daAF == 0) {
        ABBABABAcounts::AABA++;
    } else if (c.set1daAF == 1 && c.set2daAF == 1) {
        ABBABABAcounts::BBBA++;
    } else if (c.set1daAF == c.set2daAF) {
        ABBABABAcounts::p1p2++;
    } else {
        ABBABABAcounts::usedVariantsCounter++;
        // Green et al. (2010) eq. S15.2
        
        double thisDnumerator = ((1-c.set1daAF)*c.set2daAF*c.set3daAF) - (c.set1daAF*(1-c.set2daAF)*c.set3daAF);
        double thisDdenominator = ((1-c.set1daAF)*c.set2daAF*c.set3daAF) + (c.set1daAF*(1-c.set2daAF)*c.set3daAF);
        res.Dnumerator += thisDnumerator; res.lastVarsDnum += thisDnumerator; res.windowDnum += thisDnumerator;
        res.Ddenominator += thisDdenominator; res.lastVarsDdenom += thisDdenominator; res.windowDdenom += thisDdenominator;
        
        double thisF_d_denom;
        if (c.set2daAF > c.set3daAF) {
            thisF_d_denom = ((1-c.set1daAF)*c.set2daAF*c.set2daAF) - (c.set1daAF*(1-c.set2daAF)*c.set2daAF);
        } else {
            thisF_d_denom = ((1-c.set1daAF)*c.set3daAF*c.set3daAF) - (c.set1daAF*(1-c.set3daAF)*c.set3daAF);
        }
        //if (thisF_d_denom != 0) {
        /* if (thisDnumerator/thisF_d_denom > 1) {
         std::cerr << "f_d:\t" << thisDnumerator/thisF_d_denom << std::endl;
         std::cerr << "D num:\t" << thisDnumerator << std::endl;
         std::cerr << "f_d denom:\t" << thisF_d_denom << std::endl;
         std::cerr << "p1:\t" << c.set1daAF << std::endl;
         std::cerr << "p2:\t" << c.set2daAF << std::endl;
         std::cerr << "p3:\t" << c.set3daAF << std::endl;
         } */
        res.f_d_denominator += thisF_d_denom; res.window_f_d_denominator += thisF_d_denom; res.lastVarsF_d_denom += thisF_d_denom;
        res.f_d_num += thisDnumerator; res.window_f_d_num += thisDnumerator; res.lastVarsF_d_num += thisDnumerator;
        ABBABABAcounts::used_f_d_Counter++;
        //}
        
        double thisF_dM_denom;
        if (c.set1daAF <= c.set2daAF) {
            if (c.set2daAF > c.set3daAF) {
                thisF_dM_denom = ((1-c.set1daAF)*c.set2daAF*c.set2daAF) - (c.set1daAF*(1-c.set2daAF)*c.set2daAF);
            } else {
                thisF_dM_denom = ((1-c.set1daAF)*c.set3daAF*c.set3daAF) - (c.set1daAF*(1-c.set3daAF)*c.set3daAF);
            }
        } else {
            if (c.set1daAF > c.set3daAF) {
                thisF_dM_denom = -(((1-c.set1daAF)*c.set2daAF*c.set1daAF) - (c.set1daAF*(1-c.set2daAF)*c.set1daAF));
            } else {
                thisF_dM_denom = -(((1-c.set3daAF)*c.set2daAF*c.set3daAF) - (c.set3daAF*(1-c.set2daAF)*c.set3daAF));
            }
        }
        if (thisF_dM_denom == 0) {
            std::cerr << "f_d:\t" << thisDnumerator/thisF_d_denom << std::endl;
            std::cerr << "D num:\t" << thisDnumerator << std::endl;
            std::cerr << "f_d denom:\t" << thisF_d_denom << std::endl;
            std::cerr << "f_dM denom:\t" << thisF_dM_denom << std::endl;
            std::cerr << "p1:\t" << c.set1daAF << std::endl;
            std::cerr << "p2:\t" << c.set2daAF << std::endl;
            std::cerr << "p3:\t" << c.set3daAF << std::endl;
        }
        res.f_dM_denominator += thisF_dM_denom; res.window_f_dM_denominator += thisF_dM_denom; res.lastVarsF_dM_denom += thisF_dM_denom;
        
        if (c.set3daAF == 1) {
            res.f_G_denom += 1-c.set1daAF; res.f_G_num += thisDnumerator;
            res.lastVarsF_G_denom += 1-c.set1daAF; res.lastVarsF_G_num += thisDnumerator;
        }
        
        if (thisDdenominator == 0) {
            std::cerr << "P1:" << c.set1daAF << " P2:" << c.set2daAF << " P3:" << c.set3daAF << std::endl;
        }
        assert(thisDdenominator != 0);
    }
}

inline int sample01(const double& p1) {
    double r = rand() / (RAND_MAX + 1.0f);
    return r > p1;
}

inline void incrementDnumDdenomSingleSequence(ThreeSetCounts& c, ABBA_BABA_Freq_allResults& res) {
    if (c.set1daAF == -1) {
        ABBABABAcounts::noDafInfo++;
    } else if (c.set3daAF == 0) {
        ABBABABAcounts::XXAA++;
    } else {
        if (c.set3daAF == 0.5) { c.set3daAF = sample01(0.5); } // If hets then randomly sample one of the alleles
        if (c.set2daAF == 0.5) { c.set2daAF = sample01(0.5); }
        if (c.set1daAF == 0.5) { c.set1daAF = sample01(0.5); }
        
        if (c.set3daAF == 1) {
            ABBABABAcounts::XXBA++;
            // Green et al. (2010) eq. S15.1
            if (c.set1daAF == 0 && c.set2daAF == 1) {
                ABBABABAcounts::ABBA++; ABBABABAcounts::usedVariantsCounter++;
                res.Ddenominator++; res.Dnumerator++; res.lastVarsDdenom++; res.lastVarsDnum++; res.windowDdenom++; res.windowDnum++;
            } else if (c.set1daAF == 1 && c.set2daAF == 0) {
                ABBABABAcounts::BABA++; ABBABABAcounts::usedVariantsCounter++;
                res.Ddenominator++; res.Dnumerator--; res.lastVarsDdenom++; res.lastVarsDnum--; res.windowDdenom++; res.windowDnum--;
            }
        }
    }
}

inline std::string getAAfromInfo(const std::vector<std::string>& info) {
    string AA = "?";
    for (std::vector<std::string>::size_type i = 0; i != info.size(); i++) {
        string infoKey = split(info[i],'=')[0];
        if (infoKey == "AA") {
            AA = split(info[i],'=')[1];
        }
    }
    return AA;
}



void doAbbaBaba() {
    string line; // for reading the input files
    
    std::istream* vcfFile = createReader(opt::vcfFile);
    std::ifstream* setsFile = new std::ifstream(opt::setsFile.c_str());
    string setsFileRoot = stripExtension(opt::setsFile);
    std::ofstream* outFile = new std::ofstream(setsFileRoot+ "_" + opt::runName + "_abbaBaba.txt");
    string windowStartEnd = "scaffold_0\t0";
    
    // Get the sample sets
    string outgroupString; std::vector<size_t> Opos; std::vector<string> outgroup;
    if (!opt::bAaEqO) { getline(*setsFile, outgroupString); outgroup = split(outgroupString, ','); } else { outgroupString = "VCF AA field"; }
    string P3string; getline(*setsFile, P3string); std::vector<string> P3 = split(P3string, ','); std::vector<size_t> P3pos;
    string P2string; getline(*setsFile, P2string); std::vector<string> P2 = split(P2string, ','); std::vector<size_t> P2pos;
    string P1string; getline(*setsFile, P1string); std::vector<string> P1 = split(P1string, ','); std::vector<size_t> P1pos;
    if (!opt::bFrequency && (P1.size() > 1 || P2.size() > 1 || P3.size() > 1)) {
        std::cerr << "There are more than one individual on some line of the SETS.txt file" << std::endl;
        std::cerr << "Perhaps you want to use the -f option?" << std::endl;
        exit(1);
    }
    
    // Now go through the vcf and calculate D
    int totalVariantNumber = 0;
    ABBA_BABA_Freq_allResults r;
    int lastPrint = 0; int lastWindowVariant = 0;
    std::vector<double> regionDs; std::vector<double> region_f_Gs; std::vector<double> region_f_Ds; std::vector<double> region_f_DMs;
    std::vector<string> sampleNames;
    while (getline(*vcfFile, line)) {
        if (line[0] == '#' && line[1] == '#')
            continue;
        else if (line[0] == '#' && line[1] == 'C') {
            std::vector<std::string> fields = split(line, '\t');
            if (opt::sampleNameFile.empty()) {
                for (std::vector<std::string>::size_type i = NUM_NON_GENOTYPE_COLUMNS; i != fields.size(); i++) {
                    sampleNames.push_back(fields[i]);
                }
            } else {
                sampleNames = readSampleNamesFromTextFile(opt::sampleNameFile);
            }
            if (!opt::bAaEqO) { Opos = locateSet(sampleNames, outgroup); }
            P3pos = locateSet(sampleNames, P3);
            P2pos = locateSet(sampleNames, P2); P1pos = locateSet(sampleNames, P1);
            
            
            if (!opt::bAaEqO) { std::cerr << "Outgroup: "; print_vector_stream(outgroup, std::cerr); } else { std::cerr << "Outgroup: " << outgroupString << std::endl; }
            std::cerr << "P3: "; print_vector_stream(P3, std::cerr);
            std::cerr << "P2: "; print_vector_stream(P2, std::cerr);
            std::cerr << "P1: "; print_vector_stream(P1, std::cerr);
        } else {
            totalVariantNumber++;
            std::vector<std::string> fields = split(line, '\t');
            std::vector<std::string> info = split(fields[7], ';');
            if (info[0] != "INDEL") {
                if (!opt::bAaEqO) {
                    ThreeSetCounts c;
                    if (opt::bNoAaO) {
                        c = getThreeSetVariantCountsAA4(fields,P1pos,P2pos,P3pos,Opos);
                        if (opt::bFrequency) {
                            incrementDnumDdenomFrequency(c, r);
                        } else {
                            incrementDnumDdenomSingleSequence(c, r);
                        }
                    } else {
                        FourSetCounts c;
                        string AA = getAAfromInfo(info);
                        if (AA == fields[3]) {
                            c = getFourSetVariantCounts(fields,P1pos,P2pos,P3pos,Opos,"ref");
                        } else if (AA == fields[4]) {
                            c = getFourSetVariantCounts(fields,P1pos,P2pos,P3pos,Opos,"alt");
                        }
                        r.Dnumerator += ((1-c.set1daAF)*c.set2daAF*c.set3daAF*(1-c.set4daAF)) - (c.set1daAF*(1-c.set2daAF)*c.set3daAF*(1-c.set4daAF));
                        r.Ddenominator += ((1-c.set1daAF)*c.set2daAF*c.set3daAF*(1-c.set4daAF)) + (c.set1daAF*(1-c.set2daAF)*c.set3daAF*(1-c.set4daAF));
                        if (c.set2daAF > c.set3daAF) {
                            r.f_d_denominator += ((1-c.set1daAF)*c.set2daAF*c.set2daAF*(1-c.set4daAF)) - (c.set1daAF*(1-c.set2daAF)*c.set2daAF*(1-c.set4daAF));
                        } else {
                            r.f_d_denominator += ((1-c.set1daAF)*c.set3daAF*c.set3daAF*(1-c.set4daAF)) - (c.set1daAF*(1-c.set3daAF)*c.set3daAF*(1-c.set4daAF));
                        }
                    }
                } else {
                    string AA = getAAfromInfo(info);
                    ThreeSetCounts c;
                    if (AA == fields[3]) {
                        c = getThreeSetVariantCounts(fields,P1pos,P2pos,P3pos,"ref");
                    } else if (AA == fields[4]) {
                        c = getThreeSetVariantCounts(fields,P1pos,P2pos,P3pos,"alt");
                    }
                    if (opt::bFrequency) {
                        incrementDnumDdenomFrequency(c, r);
                    } else {
                        incrementDnumDdenomSingleSequence(c, r);
                    }
                }
                // if (totalVariantNumber % 100000 == 0) { std::cerr << Dnumerator << std::endl; }
            } else {
                ABBABABAcounts::indels++;
            }
            
            if (ABBABABAcounts::usedVariantsCounter % opt::windowSize == 0 && ABBABABAcounts::usedVariantsCounter != lastWindowVariant) {
                std::vector<string> s = split(windowStartEnd, '\t');
                if (s[0] == fields[0]) {
                    windowStartEnd = windowStartEnd + "\t" + fields[1];
                    if ((double)r.windowDnum/r.window_f_dM_denominator > 1) {
                        std::cerr << "D num" << r.windowDnum << std::endl;
                        std::cerr << "f_dM denom" << r.window_f_dM_denominator << std::endl;
                    }
                    *outFile << windowStartEnd << "\t" << (double)r.windowDnum/r.windowDdenom << "\t" << (double)r.window_f_d_num/r.window_f_d_denominator << "\t" << (double)r.windowDnum/r.window_f_dM_denominator << std::endl;
                    windowStartEnd = fields[0] + "\t" + fields[1];
                } else {
                    windowStartEnd = fields[0] + "\t0";
                }
                r.windowDnum = 0; r.windowDdenom = 0; r.window_f_d_num = 0; r.window_f_d_denominator = 0; r.window_f_dM_denominator = 0; lastWindowVariant = ABBABABAcounts::usedVariantsCounter;
            }
            
            
            if (ABBABABAcounts::usedVariantsCounter % opt::jackKniveWindowSize == 0 && ABBABABAcounts::usedVariantsCounter != lastPrint) {
                //if (totalVariantNumber % 100000 == 0) {
                if (opt::bFrequency)
                    assert(ABBABABAcounts::XXAA + ABBABABAcounts::AABA + ABBABABAcounts::BBBA + ABBABABAcounts::indels + ABBABABAcounts::noDafInfo + ABBABABAcounts::usedVariantsCounter + ABBABABAcounts::p1p2 == totalVariantNumber);
                if (ABBABABAcounts::usedVariantsCounter > (6 * opt::jackKniveWindowSize)) {
                    double Dstd_err = jackknive_std_err(regionDs); double f_Gstd_err = jackknive_std_err(region_f_Gs);
                    double f_Dstd_err = jackknive_std_err(region_f_Ds); double f_DMstd_err = jackknive_std_err(region_f_DMs);
                    std::cerr << totalVariantNumber << " variants processed. " << ABBABABAcounts::usedVariantsCounter << " variants used. \tD=" << (double)r.Dnumerator/r.Ddenominator << " std_err=" << Dstd_err << std::endl;
                    std::cerr << totalVariantNumber << " variants processed. " << ABBABABAcounts::usedVariantsCounter << " variants used. \tf_G=" << (double)r.f_G_num/r.f_G_denom << " std_err=" << f_Gstd_err << std::endl;
                    std::cerr << totalVariantNumber << " variants processed. " << ABBABABAcounts::usedVariantsCounter << " variants used. \tf_d=" << (double)r.f_d_num/r.f_d_denominator << " std_err=" << f_Dstd_err << std::endl;
                    std::cerr << totalVariantNumber << " variants processed. " << ABBABABAcounts::usedVariantsCounter << " variants used. \tf_dM=" << (double)r.Dnumerator/r.f_dM_denominator << " std_err=" << f_DMstd_err << std::endl;
                } else {
                    std::cerr << totalVariantNumber << " variants processed. " << ABBABABAcounts::usedVariantsCounter << " variants used. \tD=" << (double)r.Dnumerator/r.Ddenominator << std::endl;
                    std::cerr << totalVariantNumber << " variants processed. " << ABBABABAcounts::usedVariantsCounter << " variants used. \tf_G=" << (double)r.f_G_num/r.f_G_denom << std::endl;
                    std::cerr << totalVariantNumber << " variants processed. " << ABBABABAcounts::usedVariantsCounter << " variants used. \tf_d=" << (double)r.f_d_num/r.f_d_denominator << std::endl;
                    std::cerr << totalVariantNumber << " variants processed. " << ABBABABAcounts::usedVariantsCounter << " variants used. \tf_dM=" << (double)r.Dnumerator/r.f_dM_denominator << std::endl;
                }
                std::cerr << "Last used "<< opt::jackKniveWindowSize << " variants \t\t\t\tD=" << r.lastVarsDnum/r.lastVarsDdenom << std::endl;
                // std::cerr << "AAAA=" << XXAA << "; AABA=" << AABA << "; BBBA=" << BBBA << std::endl;
                regionDs.push_back(r.lastVarsDnum/r.lastVarsDdenom); region_f_Gs.push_back(r.lastVarsF_G_num/r.lastVarsF_G_denom);
                region_f_Ds.push_back(r.lastVarsF_d_num/r.lastVarsF_d_denom); region_f_DMs.push_back(r.lastVarsDnum/r.lastVarsF_dM_denom);
                r.lastVarsDnum = 0; r.lastVarsDdenom = 0; r.lastVarsF_d_num = 0; r.lastVarsF_d_denom = 0; r.lastVarsF_G_num = 0; r.lastVarsF_G_denom = 0; r.lastVarsF_dM_denom = 0;
                lastPrint = ABBABABAcounts::usedVariantsCounter;
            }
        }
    }
    
    double Dstd_err = jackknive_std_err(regionDs); double f_Gstd_err = jackknive_std_err(region_f_Gs);
    double f_Dstd_err = jackknive_std_err(region_f_Ds); double f_DMstd_err = jackknive_std_err(region_f_DMs);
    std::cerr << std::endl;
    std::cerr << totalVariantNumber << " variants processed. D=" << (double)r.Dnumerator/r.Ddenominator << " std_err=" << Dstd_err << std::endl;
    std::cerr << totalVariantNumber << " variants processed. f_G=" << (double)r.f_G_num/r.f_G_denom << " std_err=" << f_Gstd_err << std::endl;
    std::cerr << totalVariantNumber << " variants processed. f_d=" << (double)r.f_d_num/r.f_d_denominator << " std_err=" << f_Dstd_err << std::endl;
    std::cerr << totalVariantNumber << " variants processed. f_dM=" << (double)r.Dnumerator/r.f_dM_denominator << " std_err=" << f_DMstd_err << std::endl;
    
}


int abbaBabaMain(int argc, char** argv) {
    parseAbbaBabaOptions(argc, argv);
    doAbbaBaba();
    return 0;
    
}

void parseAbbaBabaOptions(int argc, char** argv) {
    bool die = false;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 's': arg >> opt::sampleNameFile; break;
            case 'f': opt::bFrequency = true; opt::jackKniveWindowSize = JACKKNIVE_WINDOW_SIZE_FREQUENCY; break;
            case 'w': arg >> opt::windowSize; break;
            case 'n': arg >> opt::runName; break;
            case OPT_AA_EQ_O: opt::bAaEqO = true; break;
            case OPT_NO_AA_O: opt::bNoAaO = true; break;
            case 'h':
                std::cout << ABBA_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    if (opt::runName == "") {
        if (opt::bFrequency) { opt::runName = "frequency"; }
        if (!opt::bFrequency) { opt::runName = "sequence"; }
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
        std::cout << "\n" << ABBA_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
}
