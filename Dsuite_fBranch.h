//
//  Dsuite_fBranch.h
//  DsuiteXcode
//
//  Created by Milan Malinsky on 11/11/2019.
//  Copyright © 2019 Milan Malinsky. All rights reserved.
//

#ifndef Dsuite_fBranch_h
#define Dsuite_fBranch_h

#include <stdio.h>
#include "Dsuite_utils.h"

void parseFbranchOptions(int argc, char** argv);
int fBranchMain(int argc, char** argv);

class Tree {
public:
    Tree(string treeString) {
        // First take care of any branch lengths
        std::regex branchLengths(":.*?(?=,|\\))");
        string treeNoBranchLengths = std::regex_replace(treeString,branchLengths,"");
        std::vector<string> tmpBranchEndNodeId;
        std::vector<string> tmpBranchStartNodeId;
        int numberOfInternalNodes = 0;
        std::regex sistersRegEx("\\(([a-zA-Z0-9_]+),([a-zA-Z0-9_]+)\\)");
        std::smatch match;
        string workingTreeCopy = treeNoBranchLengths;
        while (std::regex_search (workingTreeCopy,match,sistersRegEx)) {
            for (auto x:match) std::cout << x << " ";
            std::cout << std::endl;
            workingTreeCopy = std::regex_replace(workingTreeCopy,sistersRegEx,"internaNode#"+numToString(numberOfInternalNodes)+"X");
        }
    }
        
        
        
      /*  while working_tree_string.match(/\(([a-zA-Z0-9_]+)\,([a-zA-Z0-9_]+)\) /)
            numberOfInternalNodes += 1
            tmpBranchEndNodeId << $1
            tmpBranchAnnotation << $2
            tmpBranchStartNodeId << "internalNode#{numberOfInternalNodes}X"
            annotation_ary = $2.split(":")
            total_duration = 0
            annotation_ary.each {|i| total_duration += i.split(",")[1].to_f}
            tmpBranchDuration << total_duration
            tmpBranchEndNodeId << $3
            tmpBranchAnnotation << $4
            annotation_ary = $4.split(":")
            total_duration = 0
            annotation_ary.each {|i| total_duration += i.split(",")[1].to_f}
            tmpBranchDuration << total_duration
            tmpBranchStartNodeId << "internalNode#{numberOfInternalNodes}X"
            working_tree_string.sub!("(#{$1}:{#{$2}},#{$3}:{#{$4}})","internalNode#{numberOfInternalNodes}X")
        end
        
        */
        
        // Now process the tree
        /*std::vector<string> treeTaxonNames;
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
            } */
        
};


class Branch {
public:
    Branch(string inId, string inParentId, std::vector<string> inDaughterIds, string inTerminalSpeciesId) {
        id = inId;
        parentId = inParentId;
        if (inTerminalSpeciesId == "") {
            assert(inDaughterIds.size() == 2);
            daughterIds.push_back(inDaughterIds[0]);
            daughterIds.push_back(inDaughterIds[1]);
            terminalSpeciesId = "";
        } else {
            assert(inDaughterIds.size() == 0);
            terminalSpeciesId = inTerminalSpeciesId;
        }
    };
    
   // void getSetVariantCountsSimple(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
   // void getSetVariantCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
    
    string id;
    string parentId;
    std::vector<string> daughterIds;
    std::vector<string> progenyIds;
    string terminalSpeciesId;
    
private:
  //  void getBasicCounts(const std::vector<std::string>& genotypes, const std::map<size_t, string>& posToSpeciesMap);
};

/*

// Define a class for branches.
class Branch
    attr_reader :id, :origin, :termination, :parentId, :daughterId, :endCause, :extant, :annotation, :speciesId, :endPosition, :startPosition, :duration, :progeniesComplete, :progenyPassedOn, :progenyId, :terminalSpeciesId
    attr_writer :progeniesComplete, :progenyPassedOn, :progenyId, :terminalSpeciesId
    def initialize(id, origin, termination, parentId, daughterId, endCause, annotation)
        @id = id
        @origin = origin
        @termination = termination
        @duration = @origin - @termination
        @parentId = parentId
        @daughterId = daughterId # an array with two items
        @endCause = endCause
        @extant = false
        @annotation = annotation
        @startAnnotation = nil
        @endAnnotation = nil
        @progenyId = []
        @terminalSpeciesId = []
        @progeniesComplete = 0
        @progenyPassedOn = false
        @speciesId = "none"
    end
    def addSpeciesId(speciesId)
        @speciesId = speciesId
    end
    def updateExtant(extant)
        @extant = extant
    end
    def updateAnnotation(annotation)
        @annotation = annotation
    end
    def startAnnotation
        @annotation.split(":").last.split(",")[0]
    end
    def endAnnotation
        @annotation.split(":").first.split(",")[0]
    end
    def addStartPosition(startPosition)
        @startPosition = startPosition
    end
    def addEndPosition(endPosition)
        @endPosition = endPosition
    end
    def to_s
        string = ""
        string << "ID:                               #{@id}\n"
        string << "Origin:                           #{@origin}\n"
        string << "Termination:                      #{@termination}\n"
        string << "End cause:                        #{@endCause}\n"
        string << "Parent ID:                        #{@parentId}\n"
        string << "Daughter ID 1:                    #{@daughterId[0]}\n"
        string << "Daughter ID 2:                    #{@daughterId[1]}\n"
        string << "Species ID:                       #{@speciesId}\n"
        string << "Annotation                        #{@annotation}\n"
        string << "Position at origin:               #{@startPosition}\n"
        string << "Position at termination:          #{@endPosition}\n"
        string << "\n"
        string
    end
end

 def updateProgenyId
 # Determine the progeny of each branch (needed to know whether conditions are met, and for fossil constraints).
 # First of all, set progeniesComplete to 2 for all extinct and present branches.
 @branch.each {|b| b.progeniesComplete = 2 if b.endCause != "speciation"}
 allProgeniesComplete = false
 until allProgeniesComplete == true do
     # Set progenyPassedOn to true for the two root branches.
     @branch.each {|b| b.progenyPassedOn = true if b.parentId == "treeOrigin"}
     newlyCompleted = []
     @branch.each do |b|
         # Determine if the progeny of this branch is clear but has not been passed on to the parent yet.
         if b.progeniesComplete == 2
             newlyCompleted << b if b.progenyPassedOn == false
         end
     end
     allProgeniesComplete = true if newlyCompleted == []
     newlyCompleted.each do |b|
         # Find parent, pass progeny+self on to parents progeny, add parent.progeniesComplete += 1, and change own progenyPassedOn to true.
         @branch.each do |bb|
             if bb.id == b.parentId
                 ary = b.progenyId.dup
                 ary << b.id
                 aryFlat = ary.flatten
                 aryFlat.each {|i| bb.progenyId << i}
                 bb.progeniesComplete += 1
                 b.progenyPassedOn = true
                 break
             end
         end
     end
 end
 */
 
#endif /* Dsuite_fBranch_h */


