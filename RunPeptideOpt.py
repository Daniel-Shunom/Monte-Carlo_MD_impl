import sys
import os
import numpy as np
from math import exp, sqrt, log
import pandas as pd

from DataObjects.PeptideObject import PeptideOptData
from peptideEvaluate import peptideEvaluate, runGROMACS

from MCTree import Tree, randomlist
from Node import Node


treesearch = True
EXPLORATION_CONSTANT = 10
savetree_fname = "./tree_peptide.restart"
boltz = 8.617e-5


peptideDict = {'peptide': ['PHE', 'PHE', 'PHE']}
indata = PeptideOptData(peptideDict=peptideDict)
indata.setevalutor(peptideEvaluate)


if os.path.exists('data/current_peptide_scores.csv') is False:
    # Setting up score file
    print('Setting up score file')
    with open('data/current_peptide_scores.csv', 'w') as f:
        f.write('Peptide,init_sasa,end_sasa,AP,logP,norm_logP,norm_AP,Score\n')


def UBEnergy(nodelist, exploreconstant):
    
    # Compute the average uniqueness factor
    # Higher uniqueness, larger exploration score
    uniqscore = [1.0 for x in nodelist]
    for i, node in enumerate(nodelist):
        score = node.getuniquenessdata(nodelist=nodelist)
        uniqscore[i] = score
    uniqscore = [ item/sum(uniqscore) for item in uniqscore] 
    maxeng = max([node.getscore() for node in nodelist])


    def UCT_Unique_Score(node, uniqval, doprint=False):
        parent = node.getparent()
        energy = node.getscore()
        visits = node.getvisits()

        if parent is None:
            parenergy = node.getscore()
            parvisits = visits
        else:
            parenergy = parent.getscore()
            parvisits = parent.getvisits()
        
        depth = node.getdepth()
        playoutEList = node.getenergylist()
        childeng = [child.getscore() for child in node.getnodelist()]
        nodeEnergy = node.getscore()
        nodeweight = nodeEnergy
        avgweight = 1e300 
        cnt = 1

        if len(playoutEList) > 0:
            for energy in playoutEList:
                avgweight = min(avgweight, energy)
        avgweight = avgweight/cnt
 
        nChildren = len(node.getchildren())
        explore = 0.0
        if depth > 6:
            if doprint:
                print("Node %s (Depth %s, Visits:%s): Score:%s"%(node.getid(), depth, visits, -1e20))
            return -1e20

        if parent is None:
            explore = exploreconstant*uniqval*sqrt(log(100.0)/300.0)
        else:
            try:
                explore = exploreconstant*(uniqval*sqrt(log(parvisits)/visits))
            except (ValueError, ZeroDivisionError):
                explore =  exploreconstant*uniqval
        score = -avgweight + explore

        if doprint:
            print("Node %s (Depth %s, Visits:%s): Exploit:%s  Explore:%s Score:%s, Uniqueness:%s"%(node.getid(), depth, visits, -avgweight, explore, score, uniqval))
        return score

    keylist = {}
    for i, node in enumerate(nodelist):
        keylist[str(node)] = uniqscore[i]


    selection = sorted(nodelist, key=lambda x:x.getid())
    selection = sorted(selection, key=lambda x:UCT_Unique_Score(x, keylist[str(x)], doprint=True))[-1]
    print("Selecting Node %s with Score: %s"%(selection.getid(),  UCT_Unique_Score(selection, keylist[str(selection)], doprint=True)  ))
    return selection


nheadExpandLoop = 2 
if treesearch:

    #---Tree Class Test Code---
    print('Starting Tree Search') 
    tree = Tree(seeddata=indata, playouts=10, selectfunction=UBEnergy, headexpansion=20)
    tree.setconstant(EXPLORATION_CONSTANT)
    #tree.loadtree('tree_peptide.restart', seeddata=indata)

    for iLoop in range(1,40):
        print("Loop Number: %s"%(iLoop))

        tree.expand(nExpansions=3, writeevery=10)
        tree.simulate(nSimulations=5)
        tree.playexpand(nExpansions=5)
        tree.simulate(nSimulations=3)

        tree.savetree(savetree_fname)

        curmin = tree.getbestscore()

        if curmin < -0.1842:
            scorefile = pd.read_csv('data/current_peptide_scores.csv')
            scorefile = scorefile.sort_values(by='Score', ascending=False)
            #best_idx = scorefile.loc[scorefile['Peptide']=='SER-TYR-TYR'].index[0] + 1
            best_idx = scorefile.index[0] + 1
            best_peptide = scorefile.iloc[0]['Peptide']
            print('Found among top 3 best scoring peptide, i.e., %s with score %s after %s evaluations'%(best_peptide, -curmin, best_idx))
            print('Stopping MCTS run')
            break

        if iLoop%10==0:
            tree.setheadexpand(tree.getheadexpand()+nheadExpandLoop)
