//
//  Dsuite_utils.cpp
//  Dsuite
//
//  Created by Milan Malinsky on 02/04/2019.
//

#include "Dsuite_utils.h"

long double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
    return erfcl(-x/std::sqrt(2))/2;
}

double Fd_Denom_perVariant(double p1, double p2, double p3, double pO) {
    double Fd_Denom = 0;
    if (p2 > p3) Fd_Denom = ((1-p1)*p2*p2*(1-pO)) - (p1*(1-p2)*p2*(1-pO));
    else Fd_Denom = ((1-p1)*p3*p3*(1-pO)) - (p1*(1-p3)*p3*(1-pO));
    return Fd_Denom;
}

double fG_Denom_perVariant(double p1, double p3a, double p3b, double pO) {
    double fG_Denom = ((1-p1)*p3a*p3b*(1-pO)) - (p1*(1-p3a)*p3b*(1-pO));
    return fG_Denom;
}

double FdM_Denom_perVariant(double p1, double p2, double p3, double pO) {
    double FdM_Denom = 0;
    if (p1 <= p2) {
        if (p2 > p3) FdM_Denom = ((1-p1) * p2 * p2 * (1-pO)) - (p1 * (1-p2) * p2 * (1-pO));
        else FdM_Denom = ((1-p1) * p3 * p3 * (1-pO)) - (p1 * (1-p3) * p3 * (1-pO));
    } else {
        if (p1 > p3) FdM_Denom = -(((1-p1)*p2*p1*(1-pO)) - (p1*(1-p2)*p1*(1-pO)));
        else FdM_Denom = -(((1-p3)*p2*p3*(1-pO)) - (p3*(1-p2)*p3*(1-pO)));
    }
    return FdM_Denom;
}



// Works only on biallelic markers
void GeneralSetCounts::getSetVariantCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    
    getBasicCounts(genotypes, posToSpeciesMap);
    
    // If at least one of the outgroup individuals has non-missing data
    // Find out what is the "ancestral allele" - i.e. the one more common in the outgroup
    try {
        if (setAlleleCounts.at("Outgroup") > 0) {
            if ((double)setAltCounts.at("Outgroup")/setAlleleCounts.at("Outgroup") < 0.5) { AAint = AncestralAlleleRef; }
            else { AAint = AncestralAlleleAlt; }
        }
    } catch (std::out_of_range& e) { AAint = AncestralAlleleMissing; }
    
    // Now fill in the allele frequencies
    double totalAAF = 0; double totalDAF = 0; int numNonZeroCounts = 0;
    for(std::map<string,int>::iterator it = setAltCounts.begin(); it != setAltCounts.end(); ++it) {
        if (setAlleleCounts.at(it->first) > 0) {
            numNonZeroCounts++;
            double thisAAF = (double)setAltCounts.at(it->first)/setAlleleCounts.at(it->first);
            setAAFs[it->first] = thisAAF; totalAAF += thisAAF;
            if (AAint == 0) { // Ancestral allele seems to be the ref, so derived is alt
                setDAFs[it->first] = thisAAF; totalDAF += thisAAF;
            } else if (AAint == 1) { // Ancestral allele seems to be alt, so derived is ref
                setDAFs[it->first] = (1 - thisAAF); totalDAF += (1 - thisAAF);
            }
        }
    }
    averageAAF = totalAAF/numNonZeroCounts; averageDAF = totalDAF/numNonZeroCounts;
}

int GeneralSetCounts::returnFormatTagPosition(std::vector<std::string>& format, const std::string& tag) {
    // Find the position of GQ (genotype quality) in the genotypeData vector below
    std::vector<std::string>::iterator TAGit; int TAGi = std::numeric_limits<int>::min();
    TAGit = find (format.begin(), format.end(), tag);
    if (TAGit == format.end()) {
        // std::cerr << "This variant hasn't got associated per-sample GQ info" << std::endl;
    } else {
        TAGi = (int)std::distance( format.begin(), TAGit );
        //hasGQ = true;
    }
    return TAGi;
}


int GeneralSetCounts::checkForGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& vcfLineFields) {
    std::vector<std::string> format = split(vcfLineFields[8], ':');
    if (format.size() == 1) return LikelihoodsProbabilitiesAbsent; // The GT tag must be present in the first place
    
    int likelihoodsOrProbabilitiesTagPosition = returnFormatTagPosition(format, "GP");
    if (likelihoodsOrProbabilitiesTagPosition != std::numeric_limits<int>::min()) { likelihoodsProbabilitiesType = LikelihoodsProbabilitiesGP; }
    else {
        likelihoodsOrProbabilitiesTagPosition = returnFormatTagPosition(format, "GL");
        if (likelihoodsOrProbabilitiesTagPosition != std::numeric_limits<int>::min()) { likelihoodsProbabilitiesType = LikelihoodsProbabilitiesGL; }
        else {
            likelihoodsOrProbabilitiesTagPosition = returnFormatTagPosition(format, "PL");
            if (likelihoodsOrProbabilitiesTagPosition != std::numeric_limits<int>::min()) { likelihoodsProbabilitiesType = LikelihoodsProbabilitiesPL; }
        }
    }
    return likelihoodsOrProbabilitiesTagPosition;
}

double getExpectedGenotype(const std::vector<double>& thisProbabilities) {
    double Egenotype = thisProbabilities[1] + 2*thisProbabilities[2];
    return Egenotype;
}

void transformFromPhred(std::vector<double>& thisLikelihoods) {

    thisLikelihoods[0] = pow(10,-(thisLikelihoods[0]/10.0));
    thisLikelihoods[1] = pow(10,-(thisLikelihoods[1]/10.0));
    thisLikelihoods[2] = pow(10,-(thisLikelihoods[2]/10.0));
}

void transformFromGL(std::vector<double>& thisLikelihoods) {

    thisLikelihoods[0] = pow(10,(thisLikelihoods[0]/10.0));
    thisLikelihoods[1] = pow(10,(thisLikelihoods[1]/10.0));
    thisLikelihoods[2] = pow(10,(thisLikelihoods[2]/10.0));
}

std::vector<double> GeneralSetCounts::probabilitiesFromLikelihoods(const std::vector<double>& thisLikelihoods, const string& species) {
    std::vector<double> thisProbabilities; thisProbabilities.assign(3, 0.0);
    double multiple0 = thisLikelihoods[0]*setHWEpriorsFromAAFfromGT[species][0];
    double multiple1 = thisLikelihoods[1]*setHWEpriorsFromAAFfromGT[species][1];
    double multiple2 = thisLikelihoods[2]*setHWEpriorsFromAAFfromGT[species][2];
    double sum = multiple0 + multiple1 + multiple2;
    
    thisProbabilities[0] = multiple0/sum;
    thisProbabilities[1] = multiple1/sum;
    thisProbabilities[2] = multiple2/sum;
    
    return thisProbabilities;
}
 
void GeneralSetCounts::setHWEpriorsFromAFfromGT() {
    double AF;
    // Alternative allele frequencies
    for(std::map<string,double>::iterator it = setAAFs.begin(); it != setAAFs.end(); ++it) {
        if (it->second >= 0) AF = it->second; else AF = averageAAF; // This should be average of AFs across populations where it is known
        setHWEpriorsFromAAFfromGT[it->first][0] = pow((1-AF),2);
        setHWEpriorsFromAAFfromGT[it->first][1] = AF*(1-AF);
        setHWEpriorsFromAAFfromGT[it->first][2] = pow(AF,2);
    }
    // Derived allele frequencies
    for(std::map<string,double>::iterator it = setDAFs.begin(); it != setDAFs.end(); ++it) {
        if (it->second >= 0) AF = it->second; else AF = averageDAF; // This should be average of AFs across populations
        setHWEpriorsFromDAFfromGT[it->first][0] = pow((1-AF),2);
        setHWEpriorsFromDAFfromGT[it->first][1] = AF*(1-AF);
        setHWEpriorsFromDAFfromGT[it->first][2] = pow(AF,2);
    }
} 




void GeneralSetCounts::getAFsFromGenotypeLikelihoodsOrProbabilities(const std::vector<std::string>& genotypeFields, const std::map<size_t, string>& posToSpeciesMap, const int likelihoodsOrProbabilitiesTagPosition) {
    if (likelihoodsProbabilitiesType == LikelihoodsProbabilitiesPL || likelihoodsProbabilitiesType == LikelihoodsProbabilitiesGL) {
        setHWEpriorsFromAFfromGT();
    }
    
    for (std::vector<std::string>::size_type i = 0; i < genotypeFields.size(); i++) {
        std::string species; try { species = posToSpeciesMap.at(i); } catch (const std::out_of_range& oor) {
            continue;
        }
       // std::cerr << genotypeFields[i] << std::endl;
        std::string thisLikelihoodsOrProbabilitiesString = split(genotypeFields[i], ':')[likelihoodsOrProbabilitiesTagPosition];
        if (thisLikelihoodsOrProbabilitiesString == ".") continue;
        
        else {
            setAlleleProbCounts.at(species) += 2;
            std::vector<double> thisLikelihoodsOrProbabilities = splitToDouble(thisLikelihoodsOrProbabilitiesString,',');
            std::vector<double> thisProbabilities;
            switch (likelihoodsProbabilitiesType)
            {
                case LikelihoodsProbabilitiesPL:
                    transformFromPhred(thisLikelihoodsOrProbabilities);
                   // print_vector(thisLikelihoodsOrProbabilities, std::cerr);
                    thisProbabilities = probabilitiesFromLikelihoods(thisLikelihoodsOrProbabilities,species);
                    break;
                case LikelihoodsProbabilitiesGL: transformFromGL(thisLikelihoodsOrProbabilities);
                    thisProbabilities = probabilitiesFromLikelihoods(thisLikelihoodsOrProbabilities,species);
                    break;
                case LikelihoodsProbabilitiesGP:
                    thisProbabilities = thisLikelihoodsOrProbabilities;
                    break;
            }
            if (setAAFsFromLikelihoods.at(species) == -1) setAAFsFromLikelihoods.at(species) = 0;
            setAAFsFromLikelihoods.at(species) += getExpectedGenotype(thisProbabilities);
        }
    }
    
    for(std::map<string,double>::iterator it = setAAFsFromLikelihoods.begin(); it != setAAFsFromLikelihoods.end(); ++it) {
        if (setAAFsFromLikelihoods.at(it->first) != -1) {
            double AF = it->second/setAlleleProbCounts.at(it->first);
            it->second = AF;
            if (AAint == AncestralAlleleRef) {
                setDAFsFromLikelihoods.at(it->first) = AF;
            } else if (AAint == AncestralAlleleAlt) {
                setDAFsFromLikelihoods.at(it->first) = (1 - AF);
            }
        }
    }
     
}

void GeneralSetCountsWithSplits::getAFsFromGenotypeLikelihoodsOrProbabilitiesWithSplits(const std::vector<std::string>& genotypeFields, const std::map<size_t, string>& posToSpeciesMap, const int likelihoodsOrProbabilitiesTagPosition) {
 /* if (likelihoodsProbabilitiesType == LikelihoodsProbabilitiesPL || likelihoodsProbabilitiesType == LikelihoodsProbabilitiesGL) {
      setHWEpriorsFromAFfromGT();
  }
  
  for (std::vector<std::string>::size_type i = 0; i < genotypeFields.size(); i++) {
      std::string species; try { species = posToSpeciesMap.at(i); } catch (const std::out_of_range& oor) {
          continue;
      }
     // std::cerr << genotypeFields[i] << std::endl;
      std::string thisLikelihoodsOrProbabilitiesString = split(genotypeFields[i], ':')[likelihoodsOrProbabilitiesTagPosition];
      if (thisLikelihoodsOrProbabilitiesString == ".") continue;
      
      else {
          setAlleleProbCounts.at(species) += 2;
          std::vector<double> thisLikelihoodsOrProbabilities = splitToDouble(thisLikelihoodsOrProbabilitiesString,',');
          std::vector<double> thisProbabilities;
          switch (likelihoodsProbabilitiesType)
          {
              case LikelihoodsProbabilitiesPL:
                  transformFromPhred(thisLikelihoodsOrProbabilities);
                 // print_vector(thisLikelihoodsOrProbabilities, std::cerr);
                  thisProbabilities = probabilitiesFromLikelihoods(thisLikelihoodsOrProbabilities,species);
                  break;
              case LikelihoodsProbabilitiesGL: break;
              case LikelihoodsProbabilitiesGP:
                  thisProbabilities = thisLikelihoodsOrProbabilities;
                  break;
          }
          if (setAAFsFromLikelihoods.at(species) == -1) setAAFsFromLikelihoods.at(species) = 0;
          setAAFsFromLikelihoods.at(species) += getExpectedGenotype(thisProbabilities);
          std::cerr << abs(setAAFsFromLikelihoods.at(species)/setAlleleProbCounts.at(species)) << std::endl;
          if (abs(setAAFsFromLikelihoods.at(species)/setAlleleProbCounts.at(species)) > 1.0) {
              std::cerr << "Why is this not printing?" << std::endl;
              print_vector(thisLikelihoodsOrProbabilities, std::cerr);
              print_vector(thisProbabilities, std::cerr);
              print_vector(setHWEpriorsFromAAFfromGT.at(species), std::cerr);
              std::cerr << "AF from GT: " << setAAFs.at(species) << std::endl;
          }
      }
  }
  
  for(std::map<string,double>::iterator it = setAAFsFromLikelihoods.begin(); it != setAAFsFromLikelihoods.end(); ++it) {
      if (setAAFsFromLikelihoods.at(it->first) != -1) {
          double AF = it->second/setAlleleProbCounts.at(it->first);
          if (abs(AF) > 1.0) { std::cerr << "Something funky here: " << AF << std::endl;}
          it->second = AF;
          if (AAint == AncestralAlleleRef) {
              setDAFsFromLikelihoods.at(it->first) = AF;
          } else if (AAint == AncestralAlleleAlt) {
              setDAFsFromLikelihoods.at(it->first) = (1 - AF);
          }
      }
  } */
   
    
    if (likelihoodsProbabilitiesType == LikelihoodsProbabilitiesPL || likelihoodsProbabilitiesType == LikelihoodsProbabilitiesGL) {
        setHWEpriorsFromAFfromGT();
    }
    
    for (std::vector<std::string>::size_type i = 0; i < genotypeFields.size(); i++) {
        std::string species; try { species = posToSpeciesMap.at(i); } catch (const std::out_of_range& oor) {
            continue;
        }
       // std::cerr << genotypeFields[i] << std::endl;
        std::string thisLikelihoodsOrProbabilitiesString = split(genotypeFields[i], ':')[likelihoodsOrProbabilitiesTagPosition];
        if (thisLikelihoodsOrProbabilitiesString == ".") continue;
        
        else {
            setAlleleProbCounts.at(species) += 2;
            std::vector<double> thisLikelihoodsOrProbabilities = splitToDouble(thisLikelihoodsOrProbabilitiesString,',');
            std::vector<double> thisProbabilities;
            switch (likelihoodsProbabilitiesType)
            {
                case LikelihoodsProbabilitiesPL:
                    transformFromPhred(thisLikelihoodsOrProbabilities);
                 // print_vector(thisLikelihoodsOrProbabilities, std::cerr);
                    thisProbabilities = probabilitiesFromLikelihoods(thisLikelihoodsOrProbabilities,species);
                    break;
                case LikelihoodsProbabilitiesGL: break;
                case LikelihoodsProbabilitiesGP:
                    thisProbabilities = thisLikelihoodsOrProbabilities;
                    break;
            }
            if (setAAFsFromLikelihoods.at(species) == -1) setAAFsFromLikelihoods.at(species) = 0;
            setAAFsFromLikelihoods.at(species) += getExpectedGenotype(thisProbabilities);
     
            double r = ((double) rand() / (RAND_MAX));
            if (r < 0.5) {
                setAlleleCountsSplit1fromLikelihoods.at(species) += 2;
                if (setAAFsplit1fromLikelihoods.at(species) == -1) setAAFsplit1fromLikelihoods.at(species) = 0;
                setAAFsplit1fromLikelihoods.at(species) += getExpectedGenotype(thisProbabilities);
            } else {
                setAlleleCountsSplit2fromLikelihoods.at(species) += 2;
                if (setAAFsplit2fromLikelihoods.at(species) == -1) setAAFsplit2fromLikelihoods.at(species) = 0;
                setAAFsplit2fromLikelihoods.at(species) += getExpectedGenotype(thisProbabilities);
            }
        }
    } 
    
    for(std::map<string,double>::iterator it = setAAFsFromLikelihoods.begin(); it != setAAFsFromLikelihoods.end(); ++it) {
        if (it->second != -1) {
            double AF = it->second/setAlleleProbCounts.at(it->first);
            if (AF > 1) { std::cerr << "Something funky going on: " << it->second << "setAlleleProbCounts.at(it->first): " << setAlleleProbCounts.at(it->first)<< std::endl;}
            it->second = AF;
            if (AAint == AncestralAlleleRef) {
                setDAFsFromLikelihoods.at(it->first) = AF;
            } else if (AAint == AncestralAlleleAlt) {
                setDAFsFromLikelihoods.at(it->first) = (1 - AF);
            }
            
            if (setAAFsplit1fromLikelihoods.at(it->first) != -1) {
                double AFsplit1 = setAAFsplit1fromLikelihoods.at(it->first)/setAlleleCountsSplit1fromLikelihoods.at(it->first);
                setAAFsplit1fromLikelihoods.at(it->first) = AFsplit1;
                if (AAint == AncestralAlleleRef) {
                    setDAFsplit1fromLikelihoods.at(it->first) = AFsplit1;
                } else if (AAint == AncestralAlleleAlt) {
                    setDAFsplit1fromLikelihoods.at(it->first) = (1 - AFsplit1);
                }
            }
            
            if (setAAFsplit2fromLikelihoods.at(it->first) != -1) {
                double AFsplit2 = setAAFsplit2fromLikelihoods.at(it->first)/setAlleleCountsSplit2fromLikelihoods.at(it->first);
                setAAFsplit2fromLikelihoods.at(it->first) = AFsplit2;
                if (AAint == AncestralAlleleRef) {
                    setDAFsplit2fromLikelihoods.at(it->first) = AFsplit2;
                } else if (AAint == AncestralAlleleAlt) {
                    setDAFsplit2fromLikelihoods.at(it->first) = (1 - AFsplit2);
                }
            }
        }
    }
}



// Works only on biallelic markers
void GeneralSetCounts::getSetVariantCountsSimple(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    // std::cerr << fields[0] << "\t" << fields[1] << std::endl;
    getBasicCounts(genotypes, posToSpeciesMap);
    
    // Now fill in the allele frequencies
    for(std::map<string,int>::iterator it = setAltCounts.begin(); it != setAltCounts.end(); ++it) {
        if (setAlleleCounts.at(it->first) > 0) {
            setAAFs[it->first] = (double)setAltCounts.at(it->first)/setAlleleCounts.at(it->first);
        }
    }
}

void GeneralSetCounts::getBasicCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    // Go through the genotypes - only biallelic markers are allowed
    for (std::vector<std::string>::size_type i = 0; i != genotypes.size(); i++) {
        bool speciesDefined = true;
        std::string species; try { species = posToSpeciesMap.at(i); } catch (const std::out_of_range& oor) {
            speciesDefined = false;
        }
        // The first allele in this individual
        if (genotypes[i][0] == '1') { overall++; individualsWithVariant[i]++; }
        if (genotypes[i][2] == '1') { overall++; individualsWithVariant[i]++; }
        if (speciesDefined) {
            if (genotypes[i][0] == '1') {
                setAltCounts[species]++; setAlleleCounts[species]++;
            } else if (genotypes[i][0] == '0') {
                setAlleleCounts[species]++;
            }
            // The second allele in this individual
            if (genotypes[i][2] == '1') {
                setAltCounts[species]++; setAlleleCounts[species]++;
            } else if (genotypes[i][2] == '0') {
                setAlleleCounts[species]++;
            }
        }
    }
}

void GeneralSetCountsWithSplits::getBasicCountsWithSplits(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    // Go through the genotypes - only biallelic markers are allowed
    for (std::vector<std::string>::size_type i = 0; i != genotypes.size(); i++) {
        double r = ((double) rand() / (RAND_MAX)); bool speciesDefined = true;
        std::string species; try { species = posToSpeciesMap.at(i); } catch (const std::out_of_range& oor) {
            speciesDefined = false;
        }
        // The first allele in this individual
        if (genotypes[i][0] == '1') { overall++; individualsWithVariant[i]++; }
        if (speciesDefined) {
            if (genotypes[i][0] == '1') {
                setAltCounts[species]++; setAlleleCounts[species]++;
                if (r < 0.5) {
                    setAltCountsSplit1[species]++; setAlleleCountsSplit1[species]++;
                } else {
                    setAltCountsSplit2[species]++; setAlleleCountsSplit2[species]++;
                }
            } else if (genotypes[i][0] == '0') {
                setAlleleCounts[species]++;
                if (r < 0.5) {
                    setAlleleCountsSplit1[species]++;
                } else {
                    setAlleleCountsSplit2[species]++;
                }
            }
        }
        // The second allele in this individual
        if (genotypes[i][2] == '1') { overall++; individualsWithVariant[i]++; }
        if (speciesDefined) {
            if (genotypes[i][2] == '1') {
                setAltCounts[species]++; setAlleleCounts[species]++;
                if (r < 0.5) {
                    setAltCountsSplit1[species]++; setAlleleCountsSplit1[species]++;
                } else {
                    setAltCountsSplit2[species]++; setAlleleCountsSplit2[species]++;
                }
            } else if (genotypes[i][2] == '0') {
                setAlleleCounts[species]++;
                if (r < 0.5) {
                    setAlleleCountsSplit1[species]++;
                } else {
                    setAlleleCountsSplit2[species]++;
                }
            }
        }
    }
}

void GeneralSetCountsWithSplits::getSplitCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap) {
    
    getBasicCountsWithSplits(genotypes, posToSpeciesMap);
    
    // If at least one of the outgroup individuals has non-missing data
    // Find out what is the "ancestral allele" - i.e. the one more common in the outgroup
    try {
        if (setAlleleCounts.at("Outgroup") > 0) {
            if ((double)setAltCounts.at("Outgroup")/setAlleleCounts.at("Outgroup") < 0.5) { AAint = AncestralAlleleRef; }
            else { AAint = AncestralAlleleAlt; }
        }
    } catch (std::out_of_range& e) { AAint = -1; }
    
    // Now fill in the allele frequencies
    double totalAAF = 0; int numNonZeroCounts = 0;
    for(std::map<string,int>::iterator it = setAltCounts.begin(); it != setAltCounts.end(); ++it) {
        if (it->first == "") {
            std::cerr << "it->first " << it->first << "\t" << it->second << std::endl;
        }
        if (setAlleleCounts.at(it->first) > 0) {
            numNonZeroCounts++;
            int nSplit1; int nSplit2;
            double thisAAF = (double)setAltCounts.at(it->first)/setAlleleCounts.at(it->first);
            setAAFs[it->first] = thisAAF; totalAAF += thisAAF;
            nSplit1 = setAlleleCountsSplit1.at(it->first); nSplit2 = setAlleleCountsSplit2.at(it->first);
           // std::cerr << "it->first " << it->first << std::endl;
            try {
            if (nSplit1 > 0)
                setAAFsplit1[it->first] = (double)setAltCountsSplit1.at(it->first)/nSplit1;
            if (nSplit2 > 0)
                setAAFsplit2[it->first] = (double)setAltCountsSplit2.at(it->first)/nSplit2;
            if (AAint == AncestralAlleleRef) { // Ancestral allele seems to be the ref, so derived is alt
                setDAFs[it->first] = thisAAF;
                if (nSplit1 > 0)
                    setDAFsplit1[it->first] = (double)setAltCountsSplit1.at(it->first)/nSplit1;
                if (nSplit2 > 0)
                    setDAFsplit2[it->first] = (double)setAltCountsSplit2.at(it->first)/nSplit2;
            } else if (AAint == AncestralAlleleAlt) { // Ancestral allele seems to be alt, so derived is ref
                setDAFs[it->first] = (1 - thisAAF);
                if (nSplit1 > 0)
                    setDAFsplit1[it->first] = 1 - ((double)setAltCountsSplit1.at(it->first)/nSplit1);
                if (nSplit2 > 0)
                    setDAFsplit2[it->first] = 1 - ((double)setAltCountsSplit2.at(it->first)/nSplit2);
            }
                } catch (std::out_of_range& e) { std::cerr << "The trouble was here" << it->first << std::endl; }
        }
    }
    averageAAF = totalAAF/numNonZeroCounts;
    if (AAint == AncestralAlleleRef) averageDAF = averageAAF;
    else if (AAint == AncestralAlleleAlt) averageDAF = (1 - averageAAF);
}


double calculateOneDs(double ABBAtotal, double BABAtotal) {
    // Get the D values
    double Dnum1 = ABBAtotal - BABAtotal;
    
    double Ddenom1 = ABBAtotal + BABAtotal;
    double D = Dnum1/Ddenom1;
    return D;
}



double* calculateThreeDs(double ABBAtotal, double BABAtotal, double BBAAtotal) {
    // Get the D values
    double Dnum1 = ABBAtotal - BABAtotal;
    double Dnum2 = ABBAtotal - BBAAtotal;
    double Dnum3 = BBAAtotal - BABAtotal;
    
    double Ddenom1 = ABBAtotal + BABAtotal;
    double Ddenom2 = ABBAtotal + BBAAtotal;
    double Ddenom3 = BBAAtotal + BABAtotal;
    static double Ds[3]; Ds[0] = Dnum1/Ddenom1; Ds[1] = Dnum2/Ddenom2; Ds[2] = Dnum3/Ddenom3;
    return Ds;
}


double stringToDouble(std::string s) {
    double d;
    std::stringstream ss(s); //turn the string into a stream
    ss >> d; //convert
    return d;
}


// Remove a single file extension from the filename
std::string stripExtension(const std::string& filename)
{
    size_t suffixPos = filename.find_last_of('.');
    if(suffixPos == std::string::npos)
        return filename; // no suffix
    else
        return filename.substr(0, suffixPos);
}


void split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
    
    
void splitToDouble(const std::string &s, char delim, std::vector<double> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(stringToDouble(item));
    }
}

std::vector<double> splitToDouble(const std::string &s, char delim) {
    std::vector<double> elems;
    splitToDouble(s, delim, elems);
    return elems;
}

std::vector<std::string> split2(std::string s, string delim) {
    std::vector<std::string> elems;
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delim)) != std::string::npos) {
        token = s.substr(0, pos);
        elems.push_back(token);
        s.erase(0, pos + delim.length());
    }
    elems.push_back(s);
    return elems;
}


std::vector<size_t> locateSet(std::vector<std::string>& sample_names, const std::vector<std::string>& set) {
    std::vector<size_t> setLocs;
    for (std::vector<std::string>::size_type i = 0; i != set.size(); i++) {
        std::vector<std::string>::iterator it = std::find(sample_names.begin(), sample_names.end(), set[i]);
        if (it == sample_names.end()) {
            std::cerr << "Did not find the sample: \"" << set[i] << "\"" << std::endl;
            std::cerr << "Did not find the sample: \"" << sample_names[43] << "\"" << std::endl;
            print_vector(sample_names, std::cerr,',');
        } else {
            size_t loc = std::distance(sample_names.begin(), it);
            setLocs.push_back(loc);
        }
    }
    return setLocs;
}

//
std::string suffix(const std::string& seq, size_t len)
{
    assert(seq.length() >= len);
    return seq.substr(seq.length() - len);
}

// Returns true if the filename has an extension indicating it is compressed
bool isGzip(const std::string& filename)
{
    size_t suffix_length = sizeof(GZIP_EXT) - 1;
    
    // Assume files without an extension are not compressed
    if(filename.length() < suffix_length)
        return false;
    
    std::string extension = suffix(filename, suffix_length);
    return extension == GZIP_EXT;
}

// Ensure a filehandle is open
void assertFileOpen(std::ifstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for read\n";
        exit(EXIT_FAILURE);
    }
}
// Ensure a filehandle is open
void assertFileOpen(std::ofstream& fh, const std::string& fn)
{
    if(!fh.is_open())
    {
        std::cerr << "Error: could not open " << fn << " for write\n";
        exit(EXIT_FAILURE);
    }
}


void assertGZOpen(gzstreambase& gh, const std::string& fn)
{
    if(!gh.good())
    {
        std::cerr << "Error: could not open " << fn << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Open a file that may or may not be gzipped for reading
// The caller is responsible for freeing the handle
std::istream* createReader(const std::string& filename, std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        igzstream* pGZ = new igzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ifstream* pReader = new std::ifstream(filename.c_str(), mode);
        assertFileOpen(*pReader, filename);
        return pReader;
    }
}

// Open a file that may or may not be gzipped for writing
// The caller is responsible for freeing the handle
std::ostream* createWriter(const std::string& filename,
                           std::ios_base::openmode mode)
{
    if(isGzip(filename))
    {
        ogzstream* pGZ = new ogzstream(filename.c_str(), mode);
        assertGZOpen(*pGZ, filename);
        return pGZ;
    }
    else
    {
        std::ofstream* pWriter = new std::ofstream(filename.c_str(), mode);
        assertFileOpen(*pWriter, filename);
        return pWriter;
    }
}

bool file_exists(const std::string& name) {
    std::ifstream f(name.c_str());
    return f.good();
}


void assignSplits01FromAlleleFrequency(const double p, double& splitA, double& splitB) {
    double r = ((double) rand() / (RAND_MAX));
    if (r <= p) { splitA = 1; }
    double r2 = ((double) rand() / (RAND_MAX));
    if (r2 <= p) { splitB = 1; }
}
