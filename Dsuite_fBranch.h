//
//  Dsuite_fBranch.h
//  DsuiteXcode
//
//  Created by Milan Malinsky on 11/11/2019.
//

#ifndef Dsuite_fBranch_h
#define Dsuite_fBranch_h

#include <stdio.h>
#include "Dsuite_utils.h"
#include "Dsuite_common.h"


int fBranchMain(int argc, char** argv);
void parseFbranchOptions(int argc, char** argv);


class Branch {
public:
    Branch(string inId, string inParentId, string inDaughterId1, string inDaughterId2, string inTerminalSpeciesId) {
        id = inId;
        parentId = inParentId;
        if (inTerminalSpeciesId == "unknown") {
            daughterId1 = inDaughterId1;
            daughterId2 = inDaughterId2;
            terminalSpeciesId = "";
        } else {
            //assert(inDaughterIds.size() == 0);
            terminalSpeciesId = inTerminalSpeciesId;
            daughterId1 = "none";
            daughterId2 = "none";
        }
        progeniesComplete = 0;
        progenyPassedOn = false;
    };
    
    
    string id;
    string parentId;
    string daughterId1;
    string daughterId2;
    std::vector<string> progenyIds;
    string terminalSpeciesId;
    
    Branch* parentBranch;
    Branch* sisterBranch;
    std::vector<double> fbCvals;
    std::vector<double> ZfbCvals;
    
    int progeniesComplete;
    bool progenyPassedOn;
};


class Tree {
public:
    Tree(string treeString) {
        // First take care of any branch lengths
        std::regex branchLengths(":.*?(?=,|\\))");
        string treeNoBranchLengths = std::regex_replace(treeString,branchLengths,"");
        std::vector<string> tmpBranchEndNodeId;
        std::vector<string> tmpBranchStartNodeId;
        int numberOfInternalNodes = 0;
        std::regex sistersRegEx("\\(([a-zA-Z0-9[:s:]_-]+),([a-zA-Z0-9[:s:]_-]+)\\)");
        std::regex sistersRegExNoGroups("\\([a-zA-Z0-9[:s:]_-]+,[a-zA-Z0-9[:s:]_-]+\\)");
        std::regex comma(",");
        std::smatch match;
        string workingTreeCopy = treeNoBranchLengths;
        while (std::regex_search(workingTreeCopy,match,sistersRegEx)) {
            assert(match.size() == 3);
            // for (auto x:match) std::cout << x << " "; std::cout << std::endl;
            string nodeId = "internalNode"+numToString(numberOfInternalNodes)+"X";
            tmpBranchStartNodeId.push_back(nodeId);
            tmpBranchStartNodeId.push_back(nodeId);
            if (std::count(tmpBranchEndNodeId.begin(),tmpBranchEndNodeId.end(),match[1])) duplicateTreeValueError(match[1]);
            else tmpBranchEndNodeId.push_back(match[1]);
            if (std::count(tmpBranchEndNodeId.begin(),tmpBranchEndNodeId.end(),match[2])) duplicateTreeValueError(match[2]);
            else tmpBranchEndNodeId.push_back(match[2]);
            
            workingTreeCopy = std::regex_replace(workingTreeCopy, sistersRegExNoGroups, nodeId, std::regex_constants::format_first_only);
            // std::cout << workingTreeCopy << std::endl;
            numberOfInternalNodes++;
        }
        if (std::regex_search(workingTreeCopy,comma)) {
            std::cerr << "ERROR: The tree string could not be parsed correctly! The remaining unparsed tree string is:"  << std::endl;
            std::cerr << workingTreeCopy << std::endl;
            exit(1);
        }
        
        // Prepare arrays for temporary branch format.
        std::vector<string> tmp2BranchID;
        std::vector<string> tmp2BranchParentId;
        std::vector<string> tmp2BranchDaughterId1;
        std::vector<string> tmp2BranchDaughterId2;
        std::vector<string> tmp2BranchEndNodeId;
        
        // Prepare the first two branches in temporary format (tmpBranchEndNodeId[-1] and tmpBranchEndNodeId[-2] are the two oldest branches).
        // Test if the first root branch ends in an internal node.
        std::regex internalNodeRegEx("internalNode[0-9]+X");
        tmp2BranchID.push_back("b0");
        tmp2BranchParentId.push_back("treeOrigin");
        tmp2BranchEndNodeId.push_back(tmpBranchEndNodeId[tmpBranchEndNodeId.size()-1]);
        if (std::regex_match(tmpBranchEndNodeId[tmpBranchEndNodeId.size()-1],internalNodeRegEx)) {
            tmp2BranchDaughterId1.push_back("unborn"); tmp2BranchDaughterId2.push_back("unborn");
        } else {
            tmp2BranchDaughterId1.push_back("none"); tmp2BranchDaughterId2.push_back("none");
        }
        // Repeat the above for the second branch.
        // Test if the second root branch ends in an internal node.
        tmp2BranchID.push_back("b1");
        tmp2BranchParentId.push_back("treeOrigin");
        tmp2BranchEndNodeId.push_back(tmpBranchEndNodeId[tmpBranchEndNodeId.size()-2]);
        if (std::regex_match(tmpBranchEndNodeId[tmpBranchEndNodeId.size()-2],internalNodeRegEx)) {
            tmp2BranchDaughterId1.push_back("unborn"); tmp2BranchDaughterId2.push_back("unborn");
        } else {
            tmp2BranchDaughterId1.push_back("none"); tmp2BranchDaughterId2.push_back("none");
        }
        
        // Find out about all remaining branches until either all branches end with extinctions, or all branches have reached the present.
        int branchIdCounter = 2;
        bool treeComplete = false;
        while (!treeComplete) {
            bool change = false;
            //std::cout << "tmp2BranchID.size(): " << tmp2BranchID.size() << std::endl;
            for (int i = 0; i < tmp2BranchID.size(); i++) {
                // if a branch terminated with a speciation event in the past, then add the two daughter branches
                if (tmp2BranchDaughterId1[i] == "unborn" && tmp2BranchDaughterId2[i] == "unborn") {
                    //std::cout << "tmp2BranchEndNodeId.size(): " << tmp2BranchEndNodeId.size() << std::endl;
                    // Find the two branches that have the same start node as this branch's end node.
                    for (int j = 0; j < tmpBranchStartNodeId.size(); j++) {
                       // std::cout << "j: " << j << " i: " << i << std::endl;
                        if (tmpBranchStartNodeId[j] == tmp2BranchEndNodeId[i]) {
                            tmp2BranchID.push_back("b"+numToString(branchIdCounter));
                            //std::cout << "tmp2BranchID.size(): " << tmp2BranchID.size() << " i: " << i << std::endl;
                            tmp2BranchParentId.push_back(tmp2BranchID[i]);
                            //std::cout << "tmpBranchEndNodeId.size(): " << tmpBranchEndNodeId.size() << " j: " << j << std::endl;
                            tmp2BranchEndNodeId.push_back(tmpBranchEndNodeId[j]);
                            if (std::regex_match(tmpBranchEndNodeId[j],internalNodeRegEx)) {
                                tmp2BranchDaughterId1.push_back("unborn");
                                tmp2BranchDaughterId2.push_back("unborn");
                            } else {
                                tmp2BranchDaughterId1.push_back("none");
                                tmp2BranchDaughterId2.push_back("none");
                            }
                            // Update daughter ids of temporary parent.
                            //std::cout << "tmp2BranchDaughterId1.size(): " << tmp2BranchDaughterId1.size() << " i: " << i << std::endl;
                           // std::cout << "tmp2BranchDaughterId2.size(): " << tmp2BranchDaughterId2.size() << " i: " << i << std::endl;
                            if (tmp2BranchDaughterId1[i] == "unborn") {
                                tmp2BranchDaughterId1[i] = "b"+numToString(branchIdCounter);
                            } else {
                                tmp2BranchDaughterId2[i] = "b"+numToString(branchIdCounter);
                            }
                            // Increase the branchIdCounter
                            branchIdCounter += 1;
                            change = true;
                        }
                    }
                }
            }
            if (change == false) treeComplete = true;
        }
        
        // Fill array @branch, and at the same time, add species for terminal branches.
        std::vector<string> species;
        for (int i = 0; i < tmp2BranchID.size(); i++) {
            string speciesId;
            if (std::regex_match(tmp2BranchEndNodeId[i], internalNodeRegEx)) {
                speciesId = "unknown";
            } else {
                speciesId = tmp2BranchEndNodeId[i];
                //if (tmp2BranchParentId[i] != "treeOrigin")
                allSpecies.push_back(speciesId);
            }
            branches.push_back(new Branch(tmp2BranchID[i], tmp2BranchParentId[i], tmp2BranchDaughterId1[i], tmp2BranchDaughterId2[i], speciesId));
            
        }
    };
    
    std::vector<string> allSpecies;
    std::vector<Branch*> branches;
    void updateProgenyIds();
    void fillSisterBranches();
    
};

 


#endif /* Dsuite_fBranch_h */
