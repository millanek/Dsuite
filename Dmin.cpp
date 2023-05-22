//
//  Dmin.cpp
//  Dsuite
//
//  Created by Milan Malinsky on 02/04/2019.
//

#include "Dmin.h"
#include "Dsuite_common.h"

#define SUBPROGRAM "Dtrios"

#define DEBUG 0
#define MIN_SETS 3

static const char *DMIN_USAGE_MESSAGE =
"Usage: " PROGRAM_BIN " " SUBPROGRAM " [OPTIONS] INPUT_FILE.vcf SETS.txt\n"
"Calculate the D (ABBA/BABA) and f4-ratio statistics for all trios of species in the dataset (the outgroup being fixed)\n"
"The results are as definded in Patterson et al. 2012 (equivalent to Durand et al. 2011 when the Outgroup is fixed for the ancestral allele)\n"
"The SETS.txt should have two columns: SAMPLE_ID    SPECIES_ID\n"
"The outgroup (can be multiple samples) should be specified by using the keywork Outgroup in place of the SPECIES_ID\n"
"\n"
stdInInfo
"       -h, --help                              display this help and exit\n"
"       -k, --JKnum                             (default=20) the number of Jackknife blocks to divide the dataset into; should be at least 20 for the whole dataset\n"
"       -j, --JKwindow                          (default=NA) Jackknife block size in number of informative SNPs (as used in v0.2)\n"
"                                               when specified, this is used in place of the --JKnum option\n"
regionOption    // -r
treeOption      // -t
outOption       // -o
"       -n, --run-name                          (optional) run-name will be included in the output file name after the PREFIX\n"
"       --no-f4-ratio                           (optional) don't calculate the f4-ratio\n"
"       -l NUMLINES                             (optional) the number of lines in the VCF input - required if reading the VCF via a unix pipe\n"
"       -g, --use-genotype-probabilities        (optional) use probabilities (GP tag) or calculate them from likelihoods (GL or PL tags) using a Hardy-Weinberg prior\n"
"                                               the probabilities are used to estimate allele frequencies in each population/species\n"
"       -p, --pool-seq=MIN_DEPTH                (optional) VCF contains pool-seq data; i.e., each 'individual' is a population\n"
"                                               allele frequencies are then estimated from the AD (Allelic Depth) field, as long as there are MIN_DEPTH reads\n"
"                                               e.g MIN_DEPTH=5 may be reasonable; when there are fewer reads, the allele frequency is set to missing\n"
"       -c, --no-combine                        (optional) do not output the \"_combine.txt\" and \"_combine_stderr.txt\" files\n"
"       --KS-test-for-homoplasy                 (optional) Test whether strong ABBA-informative sites cluster along the genome\n"
"                                               use p-values output in a column called \"clustering_KS_p-val1\"\n"
//"                                               TYPE can be: 1 - clustering within a vector of all segregating sites\n"
//"                                                            2 - clustering within a vector of strong ABBA and BABA sites\n"
// "                                               TYPE=2 is less sensitive, but is robust to mutation rate variation\n"
"\n"
"\nReport bugs to " PACKAGE_BUGREPORT "\n\n";


enum { OPT_NO_F4, OPT_KS_TEST };
static const char* shortopts = "hr:n:t:j:fk:l:o:gcp:";

static const struct option longopts[] = {
    { "run-name",   required_argument, NULL, 'n' },
    { "no-combine",   required_argument, NULL, 'c' },
    { "out-prefix",   required_argument, NULL, 'o' },
    { "region",   required_argument, NULL, 'r' },
    { "tree",   required_argument, NULL, 't' },
    { "JKwindow",   required_argument, NULL, 'j' },
    { "JKnum",   required_argument, NULL, 'k' },
    { "help",   no_argument, NULL, 'h' },
    { "no-f4-ratio",   no_argument, NULL, OPT_NO_F4 },
    { "use-genotype-probabilities", no_argument, NULL, 'g'},
    { "pool-seq", required_argument, NULL, 'p'},
    { "KS-test-for-homoplasy", no_argument , NULL, OPT_KS_TEST},
    { NULL, 0, NULL, 0 }
};

namespace opt
{
    static string vcfFile;
    static string setsFile;
    static string treeFile = "";
    static string runName = "";
    static string providedOutPrefix = "";
    static int jkWindowSize = 0;
    static int jkNum = 20;
    static int regionStart = -1;
    static int regionLength = -1;
    static int providedNumLines = -1;
    static bool fStats = true;
    static bool KStest = false;
    static bool useGenotypeProbabilities = false;
    static bool poolSeq = false;
    static int poolMinDepth;
    static bool combine = true;
}


int DminMain(int argc, char** argv) {
    parseDminOptions(argc, argv);
    string line; // for reading the input files
    string outFileRoot = prepareOutFileRootString(opt::providedOutPrefix, opt::runName, opt::setsFile, opt::regionStart, opt::regionLength);
    std::istream* treeFile; std::ofstream* outFileTree;
    std::map<string,std::vector<int>> treeTaxonNamesToLoc; std::vector<int> treeLevels;
    if (opt::treeFile != "") {
        treeFile = new std::ifstream(opt::treeFile.c_str());
        if (!treeFile->good()) { std::cerr << "The file " << opt::treeFile << " could not be opened. Exiting..." << std::endl; exit(1);}
        outFileTree = new std::ofstream(outFileRoot + "_tree.txt");
        getline(*treeFile, line);
        assignTreeLevelsAndLinkToTaxa(line,treeTaxonNamesToLoc,treeLevels);
        //for (std::map<string,std::vector<int>>::iterator it = treeTaxonNamesToLoc.begin(); it != treeTaxonNamesToLoc.end(); ++it) {
        //    std::cout << "{" << it->first << "}\n";
        // }
    }
    
    int VCFlineCount = assignNumLinesToAnalyse(opt::providedNumLines, opt::regionLength, opt::vcfFile);;
    
    std::istream* vcfFile;
    if (opt::vcfFile == "stdin") {
        vcfFile = &std::cin;
    } else {
        vcfFile = createReader(opt::vcfFile.c_str());
    }
    
    // Get the sample sets
    SetInformation setInfo(opt::setsFile, MIN_SETS, OutgroupRequired);

    std::ofstream* outFileBBAA = new std::ofstream(outFileRoot+"_BBAA.txt"); assertFileOpen(*outFileBBAA, outFileRoot+"_BBAA.txt");
    std::ofstream* outFileDmin = new std::ofstream(outFileRoot+"_Dmin.txt"); assertFileOpen(*outFileDmin, outFileRoot+"_Dmin.txt");
    std::ofstream* outFileCombine; if (opt::combine) {
        outFileCombine = new std::ofstream(outFileRoot+"_combine.txt");
        assertFileOpen(*outFileCombine, outFileRoot+"_combine.txt");
    }
    std::ofstream* outFileCombineStdErr; if (opt::combine) {
        outFileCombineStdErr = new std::ofstream(outFileRoot+"_combine_stderr.txt");
        assertFileOpen(*outFileCombineStdErr, outFileRoot+"_combine_stderr.txt");
    }
    
    int nCombinations = nChoosek((int)setInfo.populations.size(),3);
    if (opt::fStats) std::cerr << "Going to calculate D and f4-ratio values for " << nCombinations << " trios" << std::endl;
    else std::cerr << "Going to calculate D values for " << nCombinations << " trios" << std::endl;
    
    if (opt::treeFile != "") { // Check that the tree contains all the populations/species
        for (int i = 0; i != setInfo.populations.size(); i++) {
            try {
                treeTaxonNamesToLoc.at(setInfo.populations[i]);
            } catch (const std::out_of_range& oor) {
                std::cerr << "Out of Range error: " << oor.what() << '\n';
                std::cerr << "setInfo.populations[i]: " << setInfo.populations[i] << '\n';
                std::cerr << CHECK_TREE_ERROR_MSG << '\n';
                exit(1);
            }
        }
    }
    
    // first, get all combinations of three sets (species):
    std::vector<std::vector<string>> trios; trios.resize(nCombinations);
    std::vector<std::vector<int>> triosInt; triosInt.resize(nCombinations);
    std::vector<bool> v(setInfo.populations.size()); std::fill(v.begin(), v.begin() + 3, true); // prepare a selection vector
    int pNum = 0;
    do {
        for (int i = 0; i < v.size(); ++i) {
            if (v[i]) { trios[pNum].push_back(setInfo.populations[i]); triosInt[pNum].push_back(i); }
        } pNum++;
    } while (std::prev_permutation(v.begin(), v.end())); // Getting all permutations of the selection vector - so it selects all combinations
    std::cerr << "Done permutations" << std::endl;
    
    // Create objects to hold the results for each trio
    std::vector<TrioDinfo> trioInfos(nCombinations); for (int i = 0; i < nCombinations; i++) { TrioDinfo info; trioInfos[i] = info; }
    
    // And need to prepare the vectors to hold allele frequency values:
    std::vector<double> allPs(setInfo.populations.size(),0.0);
    std::vector<double> allSplit1Ps(setInfo.populations.size(),0.0); std::vector<int> allSplit1Counts(setInfo.populations.size(),0);
    std::vector<double> allSplit2Ps(setInfo.populations.size(),0.0); std::vector<int> allSplit2Counts(setInfo.populations.size(),0);
    std::vector<double> allCorrectionFactors(setInfo.populations.size(),0);
    
    int totalVariantNumber = 0;
    std::vector<string> sampleNames; std::vector<std::string> fields;
    // Find out how often to report progress, based on the number of trios
    int reportProgressEvery; if (nCombinations < 1000) reportProgressEvery = 100000;
    else if (nCombinations < 100000) reportProgressEvery = 10000;
    else reportProgressEvery = 1000;
    clock_t start = clock(); clock_t startGettingCounts; clock_t startCalculation;
   // double durationGettingCounts; double durationCalculation;
    int JKblockSizeBasedOnNum = 0;
    
    //int missingLikelihoodsCount = 0;
    //int errCount = 0;
    
    while (getline(*vcfFile, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end()); // Deal with any left over \r from files prepared on Windows
        if (line[0] == '#' && line[1] == '#') {
            if (opt::regionStart == -1) { VCFlineCount--; } continue;
        } else if (line[0] == '#' && line[1] == 'C') {
            if (opt::regionStart == -1) { VCFlineCount--; } JKblockSizeBasedOnNum = (VCFlineCount/opt::jkNum)-1;
            printInitialMessageTriosQuartets(opt::regionLength, VCFlineCount, JKblockSizeBasedOnNum, opt::jkWindowSize, opt::jkNum);
            fields = split(line, '\t');
            std::vector<std::string> sampleNames(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());
            setInfo.linkSetsAndVCFpositions(sampleNames);
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
                for (int i = 0; i != trios.size(); i++) {
                    trioInfos[i].addRegionDs(P3isTrios2); trioInfos[i].addRegionDs(P3isTrios1); trioInfos[i].addRegionDs(P3isTrios0);
                }
            }
            
            if (totalVariantNumber % reportProgressEvery == 0) reportProgessVCF(totalVariantNumber, VCFlineCount, start);
            
            fields = split(line, '\t');
            std::vector<std::string> genotypes(fields.begin()+NUM_NON_GENOTYPE_COLUMNS,fields.end());

            // Only consider biallelic SNPs
            string refAllele = fields[3]; string altAllele = fields[4];
            if (refAllele.length() > 1 || altAllele.length() > 1 || altAllele == "*") {
                refAllele.clear(); refAllele.shrink_to_fit(); altAllele.clear(); altAllele.shrink_to_fit();
                genotypes.clear(); genotypes.shrink_to_fit(); continue;
            }
            
            startGettingCounts = clock();
            double p_O;
            if (opt::fStats)  {
                GeneralSetCountsWithSplits* c = new GeneralSetCountsWithSplits(setInfo.popToPosMap, (int)genotypes.size());
                c->getSplitCountsNew(genotypes, setInfo.posToPopMap);
                
                if (opt::useGenotypeProbabilities) {
                    int likelihoodsOrProbabilitiesTagPosition = c->checkForGenotypeLikelihoodsOrProbabilities(fields);
                    if (likelihoodsOrProbabilitiesTagPosition == LikelihoodsProbabilitiesAbsent) {
                        printMissingLikelihoodsWarning(fields[0], fields[1]);
                        opt::useGenotypeProbabilities = false;
                    } else c->getAFsFromGenotypeLikelihoodsOrProbabilitiesWithSplits(genotypes,setInfo.posToPopMap,likelihoodsOrProbabilitiesTagPosition, atoi(fields[1].c_str()));
                }
                
                if (opt::poolSeq) {
                    int ADtagPos = c->findADtagPosition(fields);
                    c->getAFsFromADtagWithSplits(genotypes, setInfo.popToPosMap, ADtagPos, opt::poolMinDepth);
                }
                
                p_O = c->setDAFs.at("Outgroup"); if (p_O == -1) { delete c; continue; } // We need to make sure that the outgroup is defined
                
                if (opt::useGenotypeProbabilities) {
                    for (std::vector<std::string>::size_type i = 0; i != setInfo.populations.size(); i++) {
                        try {
                            allPs[i] = c->setDAFsFromLikelihoods.at(setInfo.populations[i]);
                            allSplit1Ps[i] = c->setDAFsplit1fromLikelihoods.at(setInfo.populations[i]);
                            allSplit2Ps[i] = c->setDAFsplit2fromLikelihoods.at(setInfo.populations[i]);
                            allSplit1Counts[i] = c->setAlleleCountsSplit1fromLikelihoods.at(setInfo.populations[i]);
                            allSplit2Counts[i] = c->setAlleleCountsSplit2fromLikelihoods.at(setInfo.populations[i]);
                            if(allSplit1Ps[i] < 0) {
                                std::cerr << line << std::endl;
                            std::cerr << "setInfo.populations[i] " << setInfo.populations[i] << std::endl;
                            std::cerr << "allPs[i] " << allSplit1Ps[i] << std::endl;
                            std::cerr << "allSplit1Ps[i] " << allSplit1Ps[i] << std::endl;
                            std::cerr << "allSplit2Ps[i] " << allSplit2Ps[i] << std::endl;
                            }
                        } catch (const std::out_of_range& oor) { std::cerr << "Counts are missing some info for " << setInfo.populations[i] << std::endl; }
                    }
                   // print_vector(allPs, std::cerr);
                } else if (opt::poolSeq) {
                    for (std::vector<std::string>::size_type i = 0; i != setInfo.populations.size(); i++) {
                        try {
                                allPs[i] = c->setPoolDAFs.at(setInfo.populations[i]);
                                allSplit1Ps[i] = c->setPoolDAFsplit1.at(setInfo.populations[i]);
                                allSplit2Ps[i] = c->setPoolDAFsplit2.at(setInfo.populations[i]);
                                allSplit1Counts[i] = 1; allSplit2Counts[i] = 1;
                        } catch (const std::out_of_range& oor) { std::cerr << "Counts are missing some info for " << setInfo.populations[i] << std::endl; }
                    }
                
                } else {
                    for (std::vector<std::string>::size_type i = 0; i != setInfo.populations.size(); i++) {
                        try {
                                allPs[i] = c->setDAFs.at(setInfo.populations[i]);
                                allSplit1Ps[i] = c->setDAFsplit1.at(setInfo.populations[i]);
                                allSplit2Ps[i] = c->setDAFsplit2.at(setInfo.populations[i]);
                                allSplit1Counts[i] = c->setAlleleCountsSplit1.at(setInfo.populations[i]);
                                allSplit2Counts[i] = c->setAlleleCountsSplit2.at(setInfo.populations[i]);
                                allCorrectionFactors[i] = c->setCorrectionFactors.at(setInfo.populations[i]);
                            
                            /*if (isnan(allPs[i])) {
                                                          std::cerr << "allPs[i]: " << allPs[i] << " ; Exiting ..." << std::endl;
                                                      std::cerr << "allSplit1Ps[i]: " << allSplit1Ps[i] << " ; Exiting ..." << std::endl;
                                                      std::cerr << "allSplit2Ps[i]: " << allSplit2Ps[i] << " ; Exiting ..." << std::endl;
                                                      std::cerr << "allSplit1Counts[i]: " << allSplit1Counts[i] << " ; Exiting ..." << std::endl;
                                                      std::cerr << "allSplit2Counts[i]: " << allSplit2Counts[i] << " ; Exiting ..." << std::endl;
                                                        //  std::cerr << fields[0] << " " << fields[1] << " species[i]: " << species[i] << " ; Exiting ..." << std::endl;
                                                        //  std::cerr << genotypes[speciesToPosMap.at(species[i])[0]] << std::endl;
                                                        //  exit(1);
                                                      } */
                        } catch (const std::out_of_range& oor) { std::cerr << "Counts are missing some info for " << setInfo.populations[i] << std::endl; }
                    }
                    //print_vector(allPs, std::cerr);
                }
                delete c;
            } else {
                GeneralSetCounts* c2 = (GeneralSetCountsWithSplits*) new GeneralSetCounts(setInfo.popToPosMap, (int)genotypes.size());
                c2->getSetVariantCounts(genotypes, setInfo.posToPopMap);
                if (opt::useGenotypeProbabilities) {
                    int likelihoodsOrProbabilitiesTagPosition = c2->checkForGenotypeLikelihoodsOrProbabilities(fields);
                    if (likelihoodsOrProbabilitiesTagPosition == LikelihoodsProbabilitiesAbsent) {
                        printMissingLikelihoodsWarning(fields[0], fields[1]);
                        opt::useGenotypeProbabilities = false;
                    } else c2->getAFsFromGenotypeLikelihoodsOrProbabilities(genotypes,setInfo.posToPopMap,likelihoodsOrProbabilitiesTagPosition);
                }
                
                if (opt::poolSeq) {
                    int ADtagPos = c2->findADtagPosition(fields);
                    c2->getAFsFromADtag(genotypes,setInfo.popToPosMap,ADtagPos, opt::poolMinDepth);
                }
                
                p_O = c2->setDAFs.at("Outgroup"); if (p_O == -1) { delete c2; continue; } // We need to make sure that the outgroup is defined
                if (opt::useGenotypeProbabilities) {
                    for (std::vector<std::string>::size_type i = 0; i != setInfo.populations.size(); i++) {
                        try { allPs[i] = c2->setDAFsFromLikelihoods.at(setInfo.populations[i]); }
                        catch (const std::out_of_range& oor) { std::cerr << "Counts are missing some info for " << setInfo.populations[i] << std::endl; }
                    }
                 // print_vector(allPs, std::cerr);
                } else if (opt::poolSeq) {
                    for (std::vector<std::string>::size_type i = 0; i != setInfo.populations.size(); i++) {
                        try {allPs[i] = c2->setPoolDAFs.at(setInfo.populations[i]); }
                        catch (const std::out_of_range& oor) { std::cerr << "Counts are missing some info for " << setInfo.populations[i] << std::endl; }
                    }
                } else {
                    for (std::vector<std::string>::size_type i = 0; i != setInfo.populations.size(); i++) {
                        try {allPs[i] = c2->setDAFs.at(setInfo.populations[i]); }
                        catch (const std::out_of_range& oor) { std::cerr << "Counts are missing some info for " << setInfo.populations[i] << std::endl; }
                    }
                //print_vector(allPs, std::cerr);
                //exit(1);
                }
                delete c2;
            }
            genotypes.clear(); genotypes.shrink_to_fit();
           // durationGettingCounts = ( clock() - startGettingCounts ) / (double) CLOCKS_PER_SEC;
            
            startCalculation = clock();
            // Now calculate the D stats:
            double p_S1; double p_S2; double p_S3; double ABBA; double BABA; double BBAA; double BAAB = 0; double ABAB = 0; double AABB = 0;
            double correctionP3;
            for (int i = 0; i != trios.size(); i++) {
                p_S1 = allPs[triosInt[i][0]];
                if (p_S1 == -1) continue;  // If any member of the trio has entirely missing data, just move on to the next trio
                p_S2 = allPs[triosInt[i][1]];
                if (p_S2 == -1) continue;
                p_S3 = allPs[triosInt[i][2]];
                if (p_S3 == -1) continue;
                if (p_S1 == 0 && p_S2 == 0 && p_S3 == 0) continue; // Checking if the SNP is variable in the trio
                if (p_S1 == 1 && p_S2 == 1 && p_S3 == 1) continue; // Checking if the SNP is variable in the trio
                
                // Also no need to calculate anything if the SNP is variable in only one population
             /* if (p_S1 == 0 && p_S2 == 0 && p_O == 0) continue;
              if (p_S1 == 1 && p_S2 == 1 && p_O == 1) continue;
              if (p_S1 == 0 && p_S3 == 0 && p_O == 0) continue;
              if (p_S1 == 1 && p_S3 == 1 && p_O == 1) continue;
              if (p_S2 == 0 && p_S3 == 0 && p_O == 0) continue;
              if (p_S2 == 1 && p_S3 == 1 && p_O == 1) continue; */
                
                //std::cerr << "p_S1: " << p_S1 << " ; p_S2: " << p_S2 << " ; p_S3: " << p_S3 << std::endl;
                //std::cerr << std::endl;
                
                
                ABBA = (1-p_S1)*p_S2*p_S3*(1-p_O);
                BABA = p_S1*(1-p_S2)*p_S3*(1-p_O);
                BBAA = p_S1*p_S2*(1-p_S3)*(1-p_O);
                
                if (p_O != 0) {
                    BAAB = p_S1*(1-p_S2)*(1-p_S3)*p_O;
                    ABAB = (1-p_S1)*p_S2*(1-p_S3)*p_O;
                    AABB = (1-p_S1)*(1-p_S2)*p_S3*p_O;
                    
                    ABBA = ABBA + BAAB; BABA = BABA + ABAB; BBAA = BBAA + AABB;
                }
                
                trioInfos[i].ABBAtotal += ABBA; trioInfos[i].BABAtotal += BABA; trioInfos[i].BBAAtotal += BBAA;
                
                if (ABBA > 0.5 && (ABBA + BABA) == 0) {
                    std::cerr << "ABBA : " << ABBA << std::endl;
                    std::cerr << "BABA : " << BABA << std::endl;
                    std::cerr << "(ABBA + BABA): " << (ABBA + BABA) << std::endl;
                }
                if ((ABBA + BABA) != 0) { trioInfos[i].usedVars[0]++; trioInfos[i].totalUsedVars[0]++;
                    trioInfos[i].localD1num += ABBA - BABA; trioInfos[i].localD1denom += ABBA + BABA; }
                if ((ABBA + BBAA) != 0) { trioInfos[i].usedVars[1]++; trioInfos[i].totalUsedVars[1]++;
                    trioInfos[i].localD2num += ABBA - BBAA; trioInfos[i].localD2denom += ABBA + BBAA; }
                if ((BBAA + BABA) != 0) { trioInfos[i].usedVars[2]++; trioInfos[i].totalUsedVars[2]++;
                    trioInfos[i].localD3num += BBAA - BABA; trioInfos[i].localD3denom += BBAA + BABA; }
                
                
                if (opt::KStest) {
                    if (ABBA > 0.5) {
                       // trioInfos[i].linearStrongABBApos[0].push_back(trioInfos[i].totalUsedVars[0]);
                       // trioInfos[i].linearStrongABBApos[1].push_back(trioInfos[i].totalUsedVars[1]);
                        trioInfos[i].numStrongVars[0]++; trioInfos[i].numStrongVars[1]++;
                        trioInfos[i].linearStrongABBApos[0].push_back(totalVariantNumber);
                        trioInfos[i].linearStrongABBAposStrongSitesOnly[0].push_back(trioInfos[i].numStrongVars[0]);
                        trioInfos[i].linearStrongABBApos[1].push_back(totalVariantNumber);
                        trioInfos[i].linearStrongABBAposStrongSitesOnly[1].push_back(trioInfos[i].numStrongVars[1]);
                    }
                    if (BABA > 0.5) {
                        //trioInfos[i].linearStrongBABApos[0].push_back(trioInfos[i].totalUsedVars[0]);
                        //trioInfos[i].linearStrongBABApos[2].push_back(trioInfos[i].totalUsedVars[2]);
                        trioInfos[i].numStrongVars[0]++; trioInfos[i].numStrongVars[2]++;
                        trioInfos[i].linearStrongBABApos[0].push_back(totalVariantNumber);
                        trioInfos[i].linearStrongBABAposStrongSitesOnly[0].push_back(trioInfos[i].numStrongVars[0]);
                        trioInfos[i].linearStrongBABApos[2].push_back(totalVariantNumber);
                        trioInfos[i].linearStrongBABAposStrongSitesOnly[2].push_back(trioInfos[i].numStrongVars[2]);
                    }
                    if (BBAA > 0.5) {
                        //trioInfos[i].linearStrongABBApos[2].push_back(trioInfos[i].totalUsedVars[2]);
                        //trioInfos[i].linearStrongBABApos[1].push_back(trioInfos[i].totalUsedVars[1]);
                        trioInfos[i].numStrongVars[1]++; trioInfos[i].numStrongVars[2]++;
                        trioInfos[i].linearStrongABBApos[2].push_back(totalVariantNumber);
                        trioInfos[i].linearStrongABBAposStrongSitesOnly[2].push_back(trioInfos[i].numStrongVars[2]);
                        trioInfos[i].linearStrongBABApos[1].push_back(totalVariantNumber);
                        trioInfos[i].linearStrongBABAposStrongSitesOnly[1].push_back(trioInfos[i].numStrongVars[1]);
                    }
                }
                
                
                if (opt::fStats) {
                    
                    // f_G
                 //   int c_S1a = 0; int c_S1b = 0; int c_S2a = 0; int c_S2b = 0;int c_S3a = 0; int c_S3b = 0;
                  //  c_S3a = allSplit1Counts[triosInt[i][2]]; c_S3b = allSplit2Counts[triosInt[i][2]];
                  //  c_S2a = allSplit1Counts[triosInt[i][1]]; c_S2b = allSplit2Counts[triosInt[i][1]];
                  //  c_S1a = allSplit1Counts[triosInt[i][0]]; c_S1b = allSplit2Counts[triosInt[i][0]];
                    
                    
                    
                    double p_S1a = 0; double p_S1b = 0; double p_S2a = 0; double p_S2b = 0; double p_S3a = 0; double p_S3b = 0;
                    
                    correctionP3 = allCorrectionFactors[triosInt[i][2]];
                    
                    p_S3a = allSplit1Ps[triosInt[i][2]]; p_S3b = allSplit2Ps[triosInt[i][2]];
                    p_S2a = allSplit1Ps[triosInt[i][1]]; p_S2b = allSplit2Ps[triosInt[i][1]];
                    p_S1a = allSplit1Ps[triosInt[i][0]]; p_S1b = allSplit2Ps[triosInt[i][0]];
                    
                  //  std::cerr << "p_S1a : " << p_S1a << "; p_S1b : " << p_S1b << std::endl;
                  //  std::cerr << "p_S2a : " << p_S2a << "; p_S2b : " << p_S2b << std::endl;
                  //  std::cerr << "p_S3a : " << p_S3a << "; p_S3b : " << p_S3b << std::endl;
                    
                    assert(p_S1a >= 0); assert(p_S1b >= 0);
                    assert(p_S2a >= 0); assert(p_S2b >= 0);
                    assert(p_S3a >= 0); assert(p_S3b >= 0);
                    
                    
                    double thisFgDenom1 = fG_Denom_perVariant(p_S1,p_S3a,p_S3b,p_O);
                    double thisFgDenom1_rev = fG_Denom_perVariant(p_S2,p_S3a,p_S3b,p_O);
                    
                    trioInfos[i].F_G_denom1 += fG_Denom_perVariant(p_S1,p_S3a,p_S3b,p_O);
                    trioInfos[i].F_G_denom1_reversed += fG_Denom_perVariant(p_S2,p_S3a,p_S3b,p_O);
                    trioInfos[i].F_G_denom2 += fG_Denom_perVariant(p_S1,p_S2a,p_S2b,p_O);
                    trioInfos[i].F_G_denom2_reversed += fG_Denom_perVariant(p_S3,p_S2a,p_S2b,p_O);
                    trioInfos[i].F_G_denom3 += fG_Denom_perVariant(p_S3,p_S1a,p_S1b,p_O);
                    trioInfos[i].F_G_denom3_reversed += fG_Denom_perVariant(p_S2,p_S1a,p_S1b,p_O);
                    
                    
                    
                    
                    if (p_O != 0) {
                        thisFgDenom1 += fG_Denom_perVariant(1-p_S1,1-p_S3a,1-p_S3b,1-p_O);
                        thisFgDenom1_rev += fG_Denom_perVariant(1-p_S2,1-p_S3a,1-p_S3b,1-p_O);
                        trioInfos[i].F_G_denom1 += fG_Denom_perVariant(1-p_S1,1-p_S3a,1-p_S3b,1-p_O);
                        trioInfos[i].F_G_denom1_reversed += fG_Denom_perVariant(1-p_S2,1-p_S3a,1-p_S3b,1-p_O);
                        trioInfos[i].F_G_denom2 += fG_Denom_perVariant(1-p_S1,1-p_S2a,1-p_S2b,1-p_O);
                        trioInfos[i].F_G_denom2_reversed += fG_Denom_perVariant(1-p_S3,1-p_S2a,1-p_S2b,1-p_O);
                        trioInfos[i].F_G_denom3 += fG_Denom_perVariant(1-p_S3,1-p_S1a,1-p_S1b,1-p_O);
                        trioInfos[i].F_G_denom3_reversed += fG_Denom_perVariant(1-p_S2,1-p_S1a,1-p_S1b,1-p_O);
                    }
                    
                    /* investigating rare cases of unexpected f4-ratio values
                    if (thisFgDenom1 < 0) {
                        errCount++;
                        std::cerr << "thisFgDenom1: " << thisFgDenom1 << " ; thisFgDenom1_rev: " << thisFgDenom1_rev << std::endl;
                        std::cerr << "ABBA: " << ABBA << " ; BABA: " << BABA << " ; ABBA-BABA: " << ABBA-BABA << std::endl;
                        std::cerr << "p_S1: " << p_S1 << std::endl;
                        std::cerr << "p_S2: " << p_S2 << std::endl;
                        std::cerr << "p_S3: " << p_S3 << "; p_S3a: " << p_S3a << " ; p_S3b: " << p_S3b << std::endl;
                        std::cerr << "correctionP3: " << correctionP3 << std::endl;
                        print_vector(allPs, std::cerr);
                        print_vector(allCorrectionFactors, std::cerr);
                        std::cerr << "p_O: " << p_O << std::endl;
                        std::cerr << std::endl;
                        if (errCount > 10) {
                            exit(1);
                        }
                    
                    }
                    */
                    
                    
               /*
                // Find which topology is in agreement with the counts of the BBAA, BABA, and ABBA patterns
                                   if (BBAAtotal >= BABAtotal && BBAAtotal >= ABBAtotal) {
                                       BBAAarrangement = P3isTrios2;
                                   } else if (BABAtotal >= BBAAtotal && BABAtotal >= ABBAtotal) {
                                       BBAAarrangement = P3isTrios1;
                                   } else if (ABBAtotal >= BBAAtotal && ABBAtotal >= BABAtotal) {
                                       BBAAarrangement = P3isTrios0;
                                   }
                if (totalVariantNumber % reportProgressEvery == 0) {
                    std::cerr << trios[0][0] << "\t" << trios[0][1] << "\t" << trios[0][2] << "\n";
                    std::cerr << "p_S1a: " << p_S1a << " ; p_S1b: " << p_S1b << std::endl;
                    std::cerr << "p_S2a: " << p_S2a << " ; p_S2b: " << p_S2b << std::endl;
                    std::cerr << "p_S3a: " << p_S3a << " ; p_S3b: " << p_S3b << std::endl;
                    
                    
                    std::cerr << "ABBA-BABA: " << trioInfos[i].ABBAtotal-trioInfos[i].BABAtotal << "; ABBA - BBAA: " << trioInfos[i].ABBAtotal - trioInfos[i].BBAAtotal << "; ABBA - BBAA: " << trioInfos[i].BBAAtotal - trioInfos[i].BABAtotal << std::endl;
                    std::cerr << "trioInfos[i].F_G_denom1: " << trioInfos[i].F_G_denom1 << "; trioInfos[i].F_G_denom2: " << trioInfos[i].F_G_denom2 << "; trioInfos[i].F_G_denom3: " << trioInfos[i].F_G_denom3 << std::endl;
                    std::cerr << "trioInfos[i].F_G_denom1_reversed: " << trioInfos[i].F_G_denom1_reversed << "; trioInfos[i].F_G_denom2_reversed: " << trioInfos[i].F_G_denom2_reversed << "; trioInfos[i].F_G_denom3_reversed: " << trioInfos[i].F_G_denom3_reversed << std::endl;
                    
                    std::cerr << std::endl;
                    } */
                }
                
                // std::cerr << "trioInfos[i].localD1num" << trioInfos[i].localD1denom << std::endl;
                if (opt::jkWindowSize > 0) {
                    if (trioInfos[i].usedVars[0] == opt::jkWindowSize) { trioInfos[i].addRegionDs(P3isTrios2); }
                    if (trioInfos[i].usedVars[1] == opt::jkWindowSize) { trioInfos[i].addRegionDs(P3isTrios1); }
                    if (trioInfos[i].usedVars[2] == opt::jkWindowSize) { trioInfos[i].addRegionDs(P3isTrios0); }
                }
                // } */
            }
           // durationCalculation = ( clock() - startCalculation ) / (double) CLOCKS_PER_SEC;
        }
    }
    std::cerr << "Done processing VCF. Preparing output files..." << '\n';
    
    string header = makeHeader(false, opt::fStats, opt::KStest);
    *outFileBBAA << header << std::endl; *outFileDmin << header << std::endl;
    if (opt::treeFile != "") *outFileTree << header << std::endl;
    
    int exceptionCount = 0;
    for (int i = 0; i != trios.size(); i++) { //
        // Get the D values
        try {
            trioInfos[i].calculateFinalDs();
        } catch (const char* msg) {
            exceptionCount++;
            if (exceptionCount <= 10) {
                std::cerr << msg << std::endl;
                std::cerr << "Could not calculate p-values for the trio: " << trios[i][0] << " " << trios[i][1] << " " << trios[i][2] << std::endl;
                if (opt::jkWindowSize > 0) std::cerr << "You should probably decrease the the jackknife block size (-j option)" << std::endl;
                else std::cerr << "it looks like there aren't enough ABBA-BABA informative variants for this trio" << std::endl;
                std::cerr << std::endl;
            }
            trioInfos[i].D1_p = nan(""); trioInfos[i].D2_p = nan(""); trioInfos[i].D3_p = nan("");
        }
        
        // Find which topology is in agreement with the counts of BBAA, BABA, and ABBA
        trioInfos[i].assignBBAAarrangement();
        std::vector<string> BBAAoutVec = trioInfos[i].makeOutVec(trios[i], opt::fStats, opt::KStest, trioInfos[i].BBAAarrangement);
        print_vector(BBAAoutVec,*outFileBBAA);
        
        // Find Dmin:
        trioInfos[i].assignDminArrangement();
        std::vector<string> DminOutVec = trioInfos[i].makeOutVec(trios[i], opt::fStats, opt::KStest, trioInfos[i].DminArrangement);
        print_vector(DminOutVec,*outFileDmin);
        
        // Find which arrangement of trios is consistent with the input tree (if provided):
        if (opt::treeFile != "") {
            // Check the tree arrangement for each trio...
            if (opt::treeFile != "") {
                for (int i = 0; i != trios.size(); i++) {
                    int loc1 = treeTaxonNamesToLoc[trios[i][0]][0];
                    int loc2 = treeTaxonNamesToLoc[trios[i][1]][0];
                    int loc3 = treeTaxonNamesToLoc[trios[i][2]][0];
                    trioInfos[i].treeArrangement = trioInfos[i].assignTreeArrangement(treeLevels, loc1, loc2, loc3);
                }
            }
            std::vector<string> treeOutVec = trioInfos[i].makeOutVec(trios[i], opt::fStats, opt::KStest, trioInfos[i].treeArrangement);
            print_vector(treeOutVec,*outFileTree);
        }
        
        // Output a simple file that can be used for combining multiple local runs:
        if (opt::combine) {
            *outFileCombine << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << trioInfos[i].BBAAtotal << "\t" << trioInfos[i].BABAtotal << "\t" << trioInfos[i].ABBAtotal;
            if (opt::fStats) {
                *outFileCombine << "\t" << trioInfos[i].F_G_denom1 << "\t" << trioInfos[i].F_G_denom2 << "\t" << trioInfos[i].F_G_denom3;
                *outFileCombine << "\t" << trioInfos[i].F_G_denom1_reversed << "\t" << trioInfos[i].F_G_denom2_reversed << "\t" << trioInfos[i].F_G_denom3_reversed;
                *outFileCombine << std::endl;
            } else {
                *outFileCombine << std::endl;
            }
            print_vector(trioInfos[i].regionDs[0], *outFileCombineStdErr, ',', false); *outFileCombineStdErr << "\t"; print_vector(trioInfos[i].regionDs[1], *outFileCombineStdErr, ',', false); *outFileCombineStdErr << "\t";
            print_vector(trioInfos[i].regionDs[2], *outFileCombineStdErr, ',',false); *outFileCombineStdErr << std::endl;
        }
        //std::cerr << trios[i][0] << "\t" << trios[i][1] << "\t" << trios[i][2] << "\t" << D1 << "\t" << D2 << "\t" << D3 << "\t" << BBAAtotals[i] << "\t" << BABAtotals[i] << "\t" << ABBAtotals[i] << std::endl;
    }
    if (exceptionCount > 10) {
        std::cerr << "..." << std::endl;
        std::cerr << "p-value could not be calculated for " << exceptionCount << " trios" << std::endl;
        if (opt::jkWindowSize > 0) std::cerr << "You should probably decrease the the jackknife block size (-j option)" << std::endl;
        else std::cerr << "it looks like there aren't enough ABBA-BABA informative variants for these trios" << std::endl;
        std::cerr << "If this was a run for a subset of the genome (e.g. one chromosome), you may still get p-values for these trios from DtriosCombine" << std::endl;
        std::cerr << std::endl;
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
            case 'j': arg >> opt::jkWindowSize; break;
            case 'k': arg >> opt::jkNum; break;
            case OPT_NO_F4: opt::fStats = false; break;
            case OPT_KS_TEST: opt::KStest = true; break;
            case 'c': opt::combine = false; break;
            case 'g': opt::useGenotypeProbabilities = true; break;
            case 'l': arg >> opt::providedNumLines; break;
            case 'o': arg >> opt::providedOutPrefix; break;
            case 'p': opt::poolSeq = true; arg >> opt::poolMinDepth; break;
            case 'r': arg >> regionArgString; regionArgs = split(regionArgString, ',');
                if (regionArgs.size() != 2) {
                    std::cerr << "the --region argument should be two numbers separated by a comma\n";
                    die = true;
                } else {
                    opt::regionStart = (int)stringToDouble(regionArgs[0]); opt::regionLength = (int)stringToDouble(regionArgs[1]);  break;
                }
            case 'h':
                std::cout << DMIN_USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
        }
    }
    
    int maxNumArgs = 2; int minNumArgs = 2; // if (opt::poolSeq) { minNumArgs = 1; }
    
    if (opt::poolSeq && opt::useGenotypeProbabilities) {
        std::cerr << "Error: The -p and -g options are not compatible. Please check your command line. Exiting ....\n";
        die = true;
    }
    
    if (argc - optind < minNumArgs) {
        std::cerr << "missing arguments\n";
        die = true;
    }
    else if (argc - optind > maxNumArgs)
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
    
    if (opt::vcfFile == "stdin" && opt::providedNumLines <= 0) {
        std::cerr << "If you want to read the VCF via a pipe, you need to specify the number of lines in the input via the -l option\n";
        std::cerr << "See the example above\n";
        die = true;
    }
    
    if (die) {
        std::cout << "\n" << DMIN_USAGE_MESSAGE;
        exit(EXIT_FAILURE);
    }
}

