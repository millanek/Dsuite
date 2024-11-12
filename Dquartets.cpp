//
//  Dquartets.cpp
//  DsuiteXcode
//
//  Created by Milan Malinsky on 14/07/2020.
//

#include "Dquartets.h"
#include "Dsuite_common.h"

#define SUBPROGRAM "Dquartets"

#define DEBUG 0
#define MIN_SETS 4

static const char *DQUARTS_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf SETS.txt\n"
"Calculate the D (ABBA/BABA) and f4-ratio (f_G) statistics for all quartets of species in the dataset (there is no outgroup)\n"
"The results are as definded in Patterson et al. 2012\n"
"The SETS.txt should have two columns: SAMPLE_ID    SPECIES_ID\n"
"\n"
stdInInfo
"       -h, --help                              display this help and exit\n"
"       -k, --JKnum                             (default=20) the number of Jackknife blocks to divide the dataset into; should be at least 20 for the whole dataset\n"
"       -j, --JKwindow                          (default=NA) Jackknife block size in number of informative SNPs (as used in v0.2)\n"
"                                               when specified, this is used in place of the --JKnum option\n"
regionOption    // -r
treeOption      // -t
outOption       // -o
"       -n, --run-name                          (optional; default=quartets) run-name will be included in the output file name after the PREFIX\n"
"       --no-f4-ratio                           (optional) don't calculate the f4-ratio\n"
"       -l NUMLINES                             (optional) the number of lines in the VCF input - required if reading the VCF via a unix pipe\n"
"       -a, --allF4-ratios                      (optional) output F4 ratios for all posible arrangements\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";

enum { OPT_NO_F4 };
static const char* shortopts = "hr:n:t:j:fpk:l:o:a";

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "out-prefix",   required_argument, NULL, 'o' },
    { "region",   required_argument, NULL, 'r' },
    { "tree",   required_argument, NULL, 't' },
    { "JKwindow",   required_argument, NULL, 'j' },
    { "JKnum",   required_argument, NULL, 'k' },
    { "help",   no_argument, NULL, 'h' },
    { "no-f4-ratio",   no_argument, NULL, OPT_NO_F4 },
    { "allF4-ratios",   no_argument, NULL, 'a' },
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string treeFile = "";
    static string runName = "quartets";
    static string providedOutPrefix = "";
    static int jkWindowSize = 0;
    static int jkNum = 20;
    static int regionStart = -1;
    static int regionLength = -1;
    static int providedNumLines = -1;
    static bool fStats = true;
    static bool allF4 = false;
}


int DquartetsMain(int argc, char** argv) {
    parseDquartetsOptions(argc, argv);
    string line; // for reading the input files
    string outFileRoot = prepareOutFileRootString(opt::providedOutPrefix, opt::runName, opt::setsFile, opt::regionStart, opt::regionLength);
    
    std::istream* treeFile; std::ofstream* outFileTree;
    std::map<string,std::vector<int>> treeTaxonNamesToLoc; std::vector<int> treeLevels;
    if (opt::treeFile != "") {
        treeFile = new std::ifstream(opt::treeFile.c_str());
        if (!treeFile->good()) { std::cerr << "The file " << opt::treeFile << " could not be opened. Exiting..." << std::endl; exit(1);}
        outFileTree = new std::ofstream(outFileRoot+ "_" + opt::runName + "_tree.txt");
        getline(*treeFile, line);
        assignTreeLevelsAndLinkToTaxa(line,treeTaxonNamesToLoc,treeLevels);
        //for (std::map<string,std::vector<int>>::iterator it = treeTaxonNamesToLoc.begin(); it != treeTaxonNamesToLoc.end(); ++it) {
        //    std::cout << "{" << it->first << "}\n";
        // }
    }
    
    int VCFlineCount = assignNumLinesToAnalyse(opt::providedNumLines, opt::regionLength, opt::vcfFile);
    
    std::istream* vcfFile;
    if (opt::vcfFile == "stdin") { vcfFile = &std::cin; }
    else { vcfFile = createReader(opt::vcfFile.c_str()); }
    
    // Get the sample sets
    SetInformation setInfo(opt::setsFile, MIN_SETS, OutgroupNotRequired);
    
    std::ofstream* outFileBBAA = new std::ofstream(outFileRoot+"_BBAA.txt"); assertFileOpen(*outFileBBAA, outFileRoot+"_BBAA.txt");
    std::ofstream* outFileDmin = new std::ofstream(outFileRoot+"_Dmin.txt"); assertFileOpen(*outFileDmin, outFileRoot+"_Dmin.txt");
    std::ofstream* outFileCombine = new std::ofstream(outFileRoot+"_combine.txt"); assertFileOpen(*outFileCombine, outFileRoot+"_combine.txt");
    std::ofstream* outFileCombineStdErr = new std::ofstream(outFileRoot+"_combine_stderr.txt");
    assertFileOpen(*outFileCombineStdErr, outFileRoot+"_combine_stderr.txt");

    
    int nCombinations = nChoosek((int)setInfo.populations.size(),4);
    if (opt::fStats) std::cerr << "Going to calculate D and f4-ratio values for " << nCombinations << " quartets" << std::endl;
    else std::cerr << "Going to calculate D values for " << nCombinations << " quartets" << std::endl;
    
    if (opt::treeFile != "") { // Check that the tree contains all the populations/species
        setInfo.checkIfTreeNamesMatch(treeTaxonNamesToLoc);
    }
    
    // first, get all combinations of four sets (species):
    std::vector<std::vector<string>> quartets; quartets.resize(nCombinations);
    std::vector<std::vector<int>> quartetsInt; quartetsInt.resize(nCombinations);
    std::vector<bool> v(setInfo.populations.size()); std::fill(v.begin(), v.begin() + 4, true); // prepare a selection vector
    int pNum = 0;
    do {
        for (int i = 0; i < v.size(); ++i) {
            if (v[i]) { quartets[pNum].push_back(setInfo.populations[i]); quartetsInt[pNum].push_back(i); }
        } pNum++;
    } while (std::prev_permutation(v.begin(), v.end())); // Getting all permutations of the selection vector - so it selects all combinations
    std::cerr << "Done permutations" << std::endl;
    
    // Create objects to hold the results for each quartet
    std::vector<QuartetDinfo> quartetInfos(nCombinations); for (int i = 0; i < nCombinations; i++) {
        QuartetDinfo info; quartetInfos[i] = info;
    }
    
    // If a tree was supplied, check the tree arrangement for each trio...
    if (opt::treeFile != "") {
        for (int i = 0; i != quartets.size(); i++) {
            int loc1 = treeTaxonNamesToLoc[quartets[i][0]][0];
            int loc2 = treeTaxonNamesToLoc[quartets[i][1]][0];
            int loc3 = treeTaxonNamesToLoc[quartets[i][2]][0];
            int loc4 = treeTaxonNamesToLoc[quartets[i][3]][0];
            quartetInfos[i].treeArrangement = quartetInfos[i].assignQuartetTreeArrangement(treeLevels, loc1, loc2, loc3,loc4);
        }
    }
    
    // And need to prepare the vectors to hold allele frequency values:
    std::vector<double> allPs(setInfo.populations.size(),0.0);
    std::vector<double> allSplit1Ps(setInfo.populations.size(),0.0); std::vector<int> allSplit1Counts(setInfo.populations.size(),0);
    std::vector<double> allSplit2Ps(setInfo.populations.size(),0.0); std::vector<int> allSplit2Counts(setInfo.populations.size(),0);
    
    int totalVariantNumber = 0;
    std::vector<string> sampleNames; std::vector<std::string> fields;
    // Find out how often to report progress, based on the number of trios
    int reportProgressEvery; if (nCombinations < 1000) reportProgressEvery = 100000;
    else if (nCombinations < 100000) reportProgressEvery = 10000;
    else reportProgressEvery = 1000;
    clock_t start; clock_t startGettingCounts; clock_t startCalculation;
    double durationOverall; double durationGettingCounts; double durationCalculation;
    int JKblockSizeBasedOnNum = 0;
    
    while (getline(*vcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#') {
            VCFlineCount--; continue;
        } else if (line[0] == '#' && line[1] == 'C') {
            VCFlineCount--; JKblockSizeBasedOnNum = (VCFlineCount/opt::jkNum)-1;
            printInitialMessageTriosQuartets(opt::regionLength, VCFlineCount, JKblockSizeBasedOnNum, opt::jkWindowSize, opt::jkNum);
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            setInfo.linkSetsAndVCFpositions(sampleNames);
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
            if (totalVariantNumber % JKblockSizeBasedOnNum == 0 && opt::jkWindowSize == 0) {
                for (int i = 0; i != quartets.size(); i++) {
                    quartetInfos[i].addRegionDs(P3isTrios2); quartetInfos[i].addRegionDs(P3isTrios1); quartetInfos[i].addRegionDs(P3isTrios0);
                }
            }
            if (totalVariantNumber % reportProgressEvery == 0) {
                durationOverall = ( clock() - start ) / (double) CLOCKS_PER_SEC;
                std::cerr << "Processed " << totalVariantNumber << " variants (" << ((double)totalVariantNumber/VCFlineCount)*100 << "%) in " << durationOverall << "secs" << std::endl;
                //std::cerr << "GettingCounts " << durationGettingCounts << " calculation " << durationCalculation << "secs" << std::endl;
            }
            fields = split(line, '\t');
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());

            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4];
            if (refAllele.length() > 1 || altAllele.length() > 1 || altAllele == "*") {
                refAllele.clear(); refAllele.shrink_to_fit(); altAllele.clear(); altAllele.shrink_to_fit();
                genotypes.clear(); genotypes.shrink_to_fit(); continue;
            }
            
            startGettingCounts = clock();
            if (opt::fStats)  {
                GeneralSetCountsWithSplits* c = new GeneralSetCountsWithSplits(setInfo.popToPosMap, (int)genotypes.size());
                c->getSplitCountsNew(genotypes, setInfo.posToPopMap);
                for (std::vector<std::string>::size_type i = 0; i != setInfo.populations.size(); i++) {
                    try {
                        allPs[i] = c->setAAFs.at(setInfo.populations[i]);
                        allSplit1Ps[i] = c->setAAFsplit1.at(setInfo.populations[i]);
                        allSplit2Ps[i] = c->setAAFsplit2.at(setInfo.populations[i]);
                        allSplit1Counts[i] = c->setAlleleCountsSplit1.at(setInfo.populations[i]);
                        allSplit2Counts[i] = c->setAlleleCountsSplit2.at(setInfo.populations[i]);
                       // std::cerr << "species[i] " << species[i] << "; allPs[i] " << allPs[i] << " ; c->setDAFs[species[i]] " << c->setDAFs[0] << std::endl;
                    } catch (const std::out_of_range& oor) {
                        std::cerr << "Counts are missing some info for " << setInfo.populations[i] << std::endl;
                    }
                }
                delete c;
            } else {
                GeneralSetCounts* c = (GeneralSetCountsWithSplits*) new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size());
                c->getSetVariantCounts(genotypes, setInfo.posToPopMap);
                for (std::vector<std::string>::size_type i = 0; i != setInfo.populations.size(); i++) {
                    allPs[i] = c->setAAFs.at(setInfo.populations[i]);
                 //   std::cerr << "species[i] " << species[i] << "; allPs[i] " << allPs[i] << std::endl;
                }
                delete c;
            }
            genotypes.clear(); genotypes.shrink_to_fit();
            durationGettingCounts = ( clock() - startGettingCounts ) / (double) CLOCKS_PER_SEC;
            
            startCalculation = clock();
            // Now calculate the D stats:
            double p_S1; double p_S2; double p_S3; double p_S4; double ABBA; double BABA; double BBAA; double BAAB; double ABAB; double AABB;
            for (int i = 0; i != quartets.size(); i++) {
                p_S1 = allPs[quartetsInt[i][0]];
             //   std::cerr << "p_S1 " << p_S1 << std::endl;
                if (p_S1 == -1) continue;  // If any member of the trio has entirely missing data, just move on to the next trio
                p_S2 = allPs[quartetsInt[i][1]];
             //   std::cerr << "p_S2 " << p_S2 << std::endl;
                if (p_S2 == -1) continue;
                p_S3 = allPs[quartetsInt[i][2]];
             //   std::cerr << "p_S3 " << p_S3 << std::endl;
                if (p_S3 == -1) continue;
                p_S4 = allPs[quartetsInt[i][3]];
             //   std::cerr << "p_S4 " << p_S4 << std::endl;
                if (p_S4 == -1) continue;
                
                if (p_S1 == 0 && p_S2 == 0 && p_S3 == 0) continue; // Checking if the SNP is variable in the trio
                if (p_S1 == 0 && p_S2 == 0 && p_S4 == 0) continue; // Checking if the SNP is variable in the trio
                if (p_S1 == 0 && p_S3 == 0 && p_S4 == 0) continue; // Checking if the SNP is variable in the trio
                if (p_S2 == 0 && p_S3 == 0 && p_S4 == 0) continue; // Checking if the SNP is variable in the trio
                
                if (p_S1 == 1 && p_S2 == 1 && p_S3 == 1) continue; // Checking if the SNP is variable in the trio
                if (p_S1 == 1 && p_S2 == 1 && p_S4 == 1) continue; // Checking if the SNP is variable in the trio
                if (p_S1 == 1 && p_S3 == 1 && p_S4 == 1) continue; // Checking if the SNP is variable in the trio
                if (p_S2 == 1 && p_S3 == 1 && p_S4 == 1) continue; // Checking if the SNP is variable in the trio
                
                if (p_S4 != 1) {
                    ABBA = (1-p_S1)*p_S2*p_S3*(1-p_S4); quartetInfos[i].ABBAtotal += ABBA;
                    BABA = p_S1*(1-p_S2)*p_S3*(1-p_S4); quartetInfos[i].BABAtotal += BABA;
                    BBAA = p_S1*p_S2*(1-p_S3)*(1-p_S4); quartetInfos[i].BBAAtotal += BBAA;
                    if ((ABBA + BABA) != 0) { quartetInfos[i].usedVars[0]++; quartetInfos[i].localD1num += ABBA - BABA; quartetInfos[i].localD1denom += ABBA + BABA; }
                    if ((ABBA + BBAA) != 0) { quartetInfos[i].usedVars[1]++; quartetInfos[i].localD2num += ABBA - BBAA; quartetInfos[i].localD2denom += ABBA + BBAA; }
                    if ((BBAA + BABA) != 0) { quartetInfos[i].usedVars[2]++; quartetInfos[i].localD3num += BBAA - BABA; quartetInfos[i].localD3denom += BBAA + BABA; }
                }
                if (p_S4 != 0) {
                    BAAB = p_S1*(1-p_S2)*(1-p_S3)*p_S4; quartetInfos[i].ABBAtotal += BAAB;
                    ABAB = (1-p_S1)*p_S2*(1-p_S3)*p_S4; quartetInfos[i].BABAtotal += ABAB;
                    AABB = (1-p_S1)*(1-p_S2)*p_S3*p_S4; quartetInfos[i].BBAAtotal += AABB;
                    if (BAAB + ABAB != 0)  { quartetInfos[i].localD1num += BAAB - ABAB; quartetInfos[i].localD1denom += BAAB + ABAB; }
                    if (BAAB + AABB != 0)  { quartetInfos[i].localD2num += BAAB - AABB; quartetInfos[i].localD2denom += BAAB + AABB; }
                    if (AABB + ABAB != 0)  { quartetInfos[i].localD3num += AABB - ABAB; quartetInfos[i].localD3denom += AABB + ABAB; }
                }
                
                if (opt::fStats) {
                    
                    double p_S1a = allSplit1Ps[quartetsInt[i][0]]; double p_S1b = allSplit2Ps[quartetsInt[i][0]];
                    double p_S2a = allSplit1Ps[quartetsInt[i][1]]; double p_S2b = allSplit2Ps[quartetsInt[i][1]];
                    double p_S3a = allSplit1Ps[quartetsInt[i][2]]; double p_S3b = allSplit2Ps[quartetsInt[i][2]];
                    double p_S4a = allSplit1Ps[quartetsInt[i][3]]; double p_S4b = allSplit2Ps[quartetsInt[i][3]];
                    
                    /* Orientation 1: F4(P1, P2; P3, P4)
                     ----------------------------------
                     These are the different denominators with 'a' and 'b' being the subsamples
                     1) F4(P1, P3a; P3b, P4) ----- F_G_denom1 --- (p_S1,p_S3a,p_S3b,p_S4)
                     2) F4(P1, P2a; P2b, P4) ----- F_G_denom2 --- (p_S1,p_S2a,p_S2b,p_S4)
                     3) F4(P1a, P2; P3, P1b)
                     4) F4(P4a, P2; P3, P4b) */
                    quartetInfos[i].F_G_denoms[0] += f4_perVariant(p_S1,p_S3a,p_S3b,p_S4);
                    quartetInfos[i].F_G_denoms[1] += f4_perVariant(p_S1,p_S2a,p_S2b,p_S4);
                    quartetInfos[i].F_G_denoms[2] += f4_perVariant(p_S1a,p_S2,p_S3,p_S1b);
                    quartetInfos[i].F_G_denoms[3] += f4_perVariant(p_S4a,p_S2,p_S3,p_S4b);
                    
                    /* Orientation 1b: F4(P2, P1; P3, P4)
                    ----------------------------------   Same as Orientation 3
                    5) F4(P2, P1a; P1b, P4)
                    6) F4(P2, P3a; P3b, P4)
                    7) F4(P2a, P1; P3, P2b)
                    8) F4(P4a, P1; P3, P4b) */
                    quartetInfos[i].F_G_denoms[4] += f4_perVariant(p_S2,p_S1a,p_S1b,p_S4);
                    quartetInfos[i].F_G_denoms[5] += f4_perVariant(p_S2,p_S3a,p_S3b,p_S4);
                    quartetInfos[i].F_G_denoms[6] += f4_perVariant(p_S2a,p_S1,p_S3,p_S2b);
                    quartetInfos[i].F_G_denoms[7] += f4_perVariant(p_S4a,p_S1,p_S3,p_S4b);
                    
                    /* Orientation 2: F4(P1, P3; P2, P4)
                     ----------------------------------
                     9) F4(P1, P3a; P3b, P4) - a duplicate of 1)
                     10) F4(P1, P2a; P2b, P4) - a duplicate of 2)
                     11) F4(P1a, P3; P2, P1b)
                     12) F4(P4a, P3; P2, P4b) */
                    quartetInfos[i].F_G_denoms[8] += f4_perVariant(p_S1,p_S3a,p_S3b,p_S4);
                    quartetInfos[i].F_G_denoms[9] += f4_perVariant(p_S1,p_S2a,p_S2b,p_S4);
                    quartetInfos[i].F_G_denoms[10] += f4_perVariant(p_S1a,p_S3,p_S2,p_S1b);
                    quartetInfos[i].F_G_denoms[11] += f4_perVariant(p_S4a,p_S3,p_S2,p_S4b);
                    
                    /* Orientation 2b: F4(P3, P1; P2, P4)
                     ----------------------------------
                     13) F4(P3, P1a; P1b, P4) ---- F_G_denom3 ---   (p_S3,p_S1a,p_S1b,p_S4)
                     14) F4(P3, P2a; P2b, P4) ---- F_G_denom2_reversed --- (p_S3,p_S2a,p_S2b,p_S4)
                     15) F4(P3a, P1; P2, P3b) ----
                     16) F4(P4a, P1; P2, P4b) ---- */
                    quartetInfos[i].F_G_denoms[12] += f4_perVariant(p_S3,p_S1a,p_S1b,p_S4);
                    quartetInfos[i].F_G_denoms[13] += f4_perVariant(p_S3,p_S2a,p_S2b,p_S4);
                    quartetInfos[i].F_G_denoms[14] += f4_perVariant(p_S3a,p_S1,p_S2,p_S3b);
                    quartetInfos[i].F_G_denoms[15] += f4_perVariant(p_S4a,p_S1,p_S2,p_S4b);
                    
                    /* Orientation 3: F4(P1, P4; P2, P3)
                     ---------------------------------- Same as Orientation 1b
                     17) F4(P1, P4a; P4b, P3) - a duplicate of 8)
                     18) F4(P1a, P4; P2, P1b) - a duplicate of 5) ---- F_G_denom3_reversed --- (p_S2,p_S1a,p_S1b,p_S4)
                     19) F4(P1, P2a; P2b, P3) - a duplicate of 7)
                     20) F4(P3a, P4; P2, P3b) - a duplicate of 6) ---- F_G_denom1_reversed --- (p_S2,p_S3a,p_S3b,p_S4) */
                    quartetInfos[i].F_G_denoms[16] += f4_perVariant(p_S1,p_S4a,p_S4b,p_S3);
                    quartetInfos[i].F_G_denoms[17] += f4_perVariant(p_S1a,p_S4,p_S2,p_S1b);
                    quartetInfos[i].F_G_denoms[18] += f4_perVariant(p_S1,p_S2a,p_S2b,p_S3);
                    quartetInfos[i].F_G_denoms[19] += f4_perVariant(p_S3a,p_S4,p_S2,p_S3b);
                    
                    /* Orientation 3b: F4(P4, P1; P2, P3)
                     ----------------------------------
                     21) F4(P4, P1a; P1b, P3) - a duplicate of 13)
                     22) F4(P4, P2a; P2b, P3) - a duplicate of 14)
                     23) F4(P3a, P1; P2, P3b) - a duplicate of 15)
                     24) F4(P4a, P1; P2, P4b) - a duplicate of 16) */
                    quartetInfos[i].F_G_denoms[20] += f4_perVariant(p_S4,p_S1a,p_S1b,p_S3);
                    quartetInfos[i].F_G_denoms[21] += f4_perVariant(p_S4,p_S2a,p_S2b,p_S3);
                    quartetInfos[i].F_G_denoms[22] += f4_perVariant(p_S3a,p_S1,p_S2,p_S3b);
                    quartetInfos[i].F_G_denoms[23] += f4_perVariant(p_S4a,p_S1,p_S2,p_S4b);
                    
                   // Original version
                    quartetInfos[i].F_G_denom1 += f4_perVariant(p_S1,p_S3a,p_S3b,p_S4);
                    quartetInfos[i].F_G_denom1_reversed += f4_perVariant(p_S2,p_S3a,p_S3b,p_S4);
                    quartetInfos[i].F_G_denom2 += f4_perVariant(p_S1,p_S2a,p_S2b,p_S4);
                    quartetInfos[i].F_G_denom2_reversed += f4_perVariant(p_S3,p_S2a,p_S2b,p_S4);
                    quartetInfos[i].F_G_denom3 += f4_perVariant(p_S3,p_S1a,p_S1b,p_S4);
                    quartetInfos[i].F_G_denom3_reversed += f4_perVariant(p_S2,p_S1a,p_S1b,p_S4);
                }
                
                // std::cerr << "trioInfos[i].localD1num" << trioInfos[i].localD1denom << std::endl;
                if (opt::jkWindowSize > 0) {
                    if (quartetInfos[i].usedVars[0] == opt::jkWindowSize) { quartetInfos[i].addRegionDs(P3isTrios2); }
                    if (quartetInfos[i].usedVars[1] == opt::jkWindowSize) { quartetInfos[i].addRegionDs(P3isTrios1); }
                    if (quartetInfos[i].usedVars[2] == opt::jkWindowSize) { quartetInfos[i].addRegionDs(P3isTrios0); }
                }
                // } */
            }
            durationCalculation = ( clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
        }
    }
    std::cerr << "Done processing VCF. Preparing output files..." << '\n';
    
    string header = makeHeader(true, opt::fStats,false);
    *outFileDmin << header << std::endl;
    if(opt::allF4) {
        header += "\tF_G_denom1\tF_G_denom2\tF_G_denom3\tF_G_denom4";
    }
    if (opt::treeFile != "") *outFileTree << header << std::endl;
    *outFileBBAA << header << std::endl;
    
    int exceptionCount = 0;
    for (int i = 0; i != quartets.size(); i++) { //
        // Get the D values
        try {
            /*std::cerr << "Here..." << '\n';
            std::cerr << "quartetInfos[i]." << quartetInfos[i].ABBAtotal << '\n';
            std::cerr << "quartetInfos[i]." << quartetInfos[i].BBAAtotal << '\n';
            std::cerr << "quartetInfos[i]." << quartetInfos[i].BABAtotal << '\n'; */
            quartetInfos[i].calculateFinalDs();
        } catch (const char* msg) {
            exceptionCount++;
            if (exceptionCount <= 10) {
                std::cerr << msg << std::endl;
                std::cerr << "Could not calculate p-values for the quartet: " << quartets[i][0] << " " << quartets[i][1] << " " << quartets[i][2] << " " << quartets[i][3]<< std::endl;
                if (opt::jkWindowSize > 0) std::cerr << "You should probably decrease the the jackknife block size (-j option)" << std::endl;
                else std::cerr << "it looks like there aren't enough ABBA-BABA informative variants for this quartet" << std::endl;
                std::cerr << std::endl;
            }
            quartetInfos[i].D1_p = nan(""); quartetInfos[i].D2_p = nan(""); quartetInfos[i].D3_p = nan("");
        }
       // std::cerr << "Here..." << '\n';
        
        // Find which topology is in agreement with the counts of BBAA, BABA, and ABBA
        quartetInfos[i].assignBBAAarrangement();
        std::vector<string> BBAAoutVec = quartetInfos[i].makeOutVec(quartets[i], opt::fStats, quartetInfos[i].BBAAarrangement, opt::allF4);
       // std::cerr << "quartetInfos[i].BBAAarrangement: " << quartetInfos[i].BBAAarrangement << std::endl;
        print_vector(BBAAoutVec,*outFileBBAA);
        
        // Find Dmin:
        quartetInfos[i].assignDminArrangement();
       // std::cerr << "quartetInfos[i].DminArrangement " << quartetInfos[i].DminArrangement << '\n';
        std::vector<string> DminOutVec = quartetInfos[i].makeOutVec(quartets[i], opt::fStats, quartetInfos[i].DminArrangement);
        print_vector(DminOutVec,*outFileDmin);
        
        // Find which arrangement of trios is consistent with the input tree (if provided):
        if (opt::treeFile != "") {
       //     std::cerr << "quartetInfos[i].treeArrangement " << quartetInfos[i].treeArrangement << '\n';
            std::vector<string> treeOutVec = quartetInfos[i].makeOutVec(quartets[i], opt::fStats, quartetInfos[i].treeArrangement, opt::allF4);
            print_vector(treeOutVec,*outFileTree);
        }
        
        // Output a simple file that can be used for combining multiple local runs:
        *outFileCombine << quartets[i][0] << "\t" << quartets[i][1] << "\t" << quartets[i][2] << "\t" << quartetInfos[i].BBAAtotal << "\t" << quartetInfos[i].BABAtotal << "\t" << quartetInfos[i].ABBAtotal;
        if (opt::fStats) {
            *outFileCombine << "\t" << quartetInfos[i].F_G_denom1 << "\t" << quartetInfos[i].F_G_denom2 << "\t" << quartetInfos[i].F_G_denom3;
            *outFileCombine << "\t" << quartetInfos[i].F_G_denom1_reversed << "\t" << quartetInfos[i].F_G_denom2_reversed << "\t" << quartetInfos[i].F_G_denom3_reversed;
            *outFileCombine << std::endl;
        } else {
            *outFileCombine << std::endl;
        }
        print_vector(quartetInfos[i].regionDs[0], *outFileCombineStdErr, ',', false); *outFileCombineStdErr << "\t"; print_vector(quartetInfos[i].regionDs[1], *outFileCombineStdErr, ',', false); *outFileCombineStdErr << "\t";
        print_vector(quartetInfos[i].regionDs[2], *outFileCombineStdErr, ',',false); *outFileCombineStdErr << std::endl;
        
        //std::cerr << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << D1 << "\t" << D2 << "\t" << D3 << "\t" << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
    }
    if (exceptionCount > 10) {
        std::cerr << "..." << std::endl;
        std::cerr << "p-value could not be calculated for " << exceptionCount << " quartets" << std::endl;
        if (opt::jkWindowSize > 0) std::cerr << "You should probably decrease the the jackknife block size (-j option)" << std::endl;
        else std::cerr << "it looks like there aren't enough ABBA-BABA informative variants for these quartets" << std::endl;
       // std::cerr << "If this was a run for a subset of the genome (e.g. one chromosome), you may still get p-values for these quartets from DtriosCombine" << std::endl;
        std::cerr << std::endl;
    }
    return 0;
    
}


void parseDquartetsOptions(int argc, char** argv) {
    bool die = false; string regionArgString; std::vector<string> regionArgs;
    for (char c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;)
    {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c)
        {
            case '?': die = true; break;
            case 'n': arg >> opt::runName; break;
            case 't': arg >> opt::treeFile; break;
            case 'j': arg >> opt::jkWindowSize; break;
            case 'k': arg >> opt::jkNum; break;
            case 'o': arg >> opt::providedOutPrefix; break;
            case OPT_NO_F4: opt::fStats = false; break;
            case 'l': arg >> opt::providedNumLines; break;
            case 'a': opt::allF4 = true; break;
            case 'r': arg >> regionArgString; regionArgs = split(regionArgString, ',');
                opt::regionStart = (int)stringToDouble(regionArgs[0]); opt::regionLength = (int)stringToDouble(regionArgs[1]);  break;
            case 'h':
                std::cout << DQUARTS_USAGE_MESSAGE;
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
        std::cout << "\n" << DQUARTS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
    
    // Parse the input filenames
    opt::vcfFile = argv[optind++];
    opt::setsFile = argv[optind++];
    
    if (opt::vcfFile == "stdin" && opt::providedNumLines <= 0) {
        std::cerr << "If you want to read the VCF via a pipe, you need to specify the number of lines in the input via the -l option\n";
        std::cerr << "See the example above\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DQUARTS_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}
