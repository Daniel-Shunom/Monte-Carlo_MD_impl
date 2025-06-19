import random
import numpy as np 
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit import DataStructs
import copy
import pandas as pd

from peptideEvaluate import que_playout_evaluate, que_playout_results
from peptideUtilities import get_short_sequence

import sys

#------------------------------------
# User Variable to Change

que_submit = False
seqRange=[1,21]
peptideNames = {'1': 'ARG',
                '2': 'HIS',
                '3': 'LYS',
                '4': 'ASP',
                '5': 'GLU',
                '6': 'SER',
                '7': 'THR',
                '8': 'ASN',
                '9': 'GLN',
                '10': 'CYS',
                '11': 'PHE', #DUPLICATE
                '12': 'GLY',
                '13': 'PRO',
                '14': 'ALA',
                '15': 'ILE',
                '16': 'LEU',
                '17': 'MET',
                '18': 'PHE', #DUPLICATE
                '19': 'TRP',
                '20': 'TYR',
                '21': 'VAL'}

peptide_shortSeq = {'ARG': 'R',
                'HIS': 'H',
                'LYS': 'K',
                'ASP': 'D',
                'GLU': 'E',
                'SER': 'S',
                'THR': 'T',
                'ASN': 'N',
                'GLN': 'Q',
                'CYS': 'C',
                'GLY': 'G',
                'PRO': 'P',
                'ALA': 'A',
                'ILE': 'I',
                'LEU': 'L',
                'MET': 'M',
                'PHE': 'F',
                'TRP': 'W',
                'TYR': 'Y',
                'VAL': 'V'}

#pepStructuredir='/home/rbatra/Research/peptide_sequence/5mer/2_Creating_coordinates/'
#pepStructuredir='/Users/rbatra/Research/Research_ANL/covid19/srilok/peptide/pdb/'
#pepStructuredir='/home/rbatra/peptide_assembly/5mer/2_Creating_coordinates/'

n_ML_playouts = 5
#------------------------------------------------------------


'''
.. module::  TemplateObject
   :platform: Unix, Windows, Linux
   :synopsis:  Data Object "Skeleton" with the functions that must be defined to properly run within the MCTS frame work. 

.. moduleauthor:: Troy D Loeffler <tloeffler@anl.gov>
.. note::  Data Object "Skeleton" with the functions that must be defined to properly run within the MCTS frame work. This outlines the basic functions
 required for a Data Object to work properly with the upper level MCTS code.

 This object's purpose is to define the opperations that are unique to a given problem.  For example a parameter search is going to require
 different opperations than optimization on a discrete grid. 

 Terminology:
   Data Object => Python Class which acts as a go between wrapper from the high level MCTS calls to the low level problem specific code.

   Data Set => Problem specific variables contained within the Data Object. For example in a parameter optimization this is the list of problem parameters (Xi = [x1, x2, x3...xn]) that are being optimized for. 

   Objective Score => Score returned by the user defined objective function that the MCTs is attempting to Minimize/Maximize (IE a f(Xi) function)                      Also known as a cost function, score function, target function, etc.
'''


class PeptideOptData(object):
    '''
     Example Data Object Class which shows the functions that must be defined
     for an object to be compatible with the MCTree
    '''
    #------------------------------------------------

    def __init__(self, peptideDict=None, 
            seqRange=[1,21],                 
            peptideNames = {'1': 'ARG',
                            '2': 'HIS',
                            '3': 'LYS',
                            '4': 'ASP',
                            '5': 'GLU',
                            '6': 'SER',
                            '7': 'THR',
                            '8': 'ASN',
                            '9': 'GLN',
                            '10': 'CYS',
                            '11': 'PHE', #DUPLICATE
                            '12': 'GLY',
                            '13': 'PRO',
                            '14': 'ALA',
                            '15': 'ILE',
                            '16': 'LEU',
                            '17': 'MET',
                            '18': 'PHE', #DUPLICATE
                            '19': 'TRP',
                            '20': 'TYR',
                            '21': 'VAL'},                 
            peptide=None, 
            evaluator=None):
        
        if peptideDict: 
            self.peptide = peptideDict['peptide']
            self.nblocks = len(self.peptide)
            self.fingerprint = None
            self.smiles = None
            self.evaluator = None
            self.short_sequence = get_short_sequence(self.peptide)
        else: 
            self.peptide = None
            self.nblocks = 3 
            self.fingerprint = None
            self.evaluator = None
            self.smiles = None            
            self.makeRandomPeptide(peptideNames, seqRange)
            self.short_sequence = get_short_sequence(self.peptide)
            

    def _sampleGivenRange(self, ParamRange): 
        lb = ParamRange[0]
        ub = ParamRange[1]            
        return random.randint(lb,ub)

    
    def makeRandomPeptide(self, peptideNames, seqRange):                 
        print('Making Random Peptide')
        #SatisfiesConstrains = False

        '''
        ##### Old code ##########
        self.a = peptideNames[str(self._sampleGivenRange(seqRange))]
        self.b = peptideNames[str(self._sampleGivenRange(seqRange))]
        self.c = peptideNames[str(self._sampleGivenRange(seqRange))]
        self.peptide = [self.a, self.b, self.c]
        self.nblocks = len(self.peptide)
        #########################
        '''

        for i in range(self.nblocks):
            self.peptide.append(peptideNames[str(self._sampleGivenRange(seqRange))])

        
        
    def setevalutor(self,evalFunc): 
        self.evaluator = evalFunc


        
    def getPeptideDict(self): 
        peptideDict={} 
        peptideDict['peptide'] = self.peptide
        peptideDict['nblocks'] = len(self.peptide)
        return peptideDict


    
    def mutate_sequence(self, nodeDepth=None):
        peptideDict = copy.deepcopy(self.getPeptideDict().copy())
        nblocks = peptideDict['nblocks']
        #newpeptide = peptideDict['peptide']
        #mutatePos = random.randint(1,len(newpeptide)) - 1
        #mutateString = peptideNames[str(self._sampleGivenRange(seqRange))]
        #newpeptide[mutatePos] = mutateString

        if nodeDepth == 0: # Random Sequence for all blocks
            newpeptide = []
            mutatePos = 'All'
            mutateString = 'XXX'
            for i in range(nblocks):
                newpeptide.append(peptideNames[str(self._sampleGivenRange(seqRange))])


        if nodeDepth > 0:  # Random Sequence for  subset blocks
            newpeptide = peptideDict['peptide']
            mutateString = None
            mutatePos = random.sample(range(nblocks), max(nblocks - nodeDepth, 1))
            for pos in mutatePos:
                mutateString = peptideNames[str(self._sampleGivenRange(seqRange))]
                newpeptide[pos] = mutateString
        
        #print('Old Peptide, New mutatePos, mutateString, Peptide: ', self.getPeptideDict()['peptide'], mutatePos, mutateString, newpeptide)
        peptideDict['peptide'] = newpeptide
        
        return peptideDict    

        
        
    def perturbate(self, node=None, inset=None, setOnly=False):
        '''
        Perturbate is the function called by the MCTS to generate a new node or playout structure.  
        The purpose of this function is to take the information within this data object and mutate it
        slightly.  This can then return a new data object of the same type

        In Signature 
          node => The node who owns this data object is passed in so that the data object may obtain information related to it's tree positioning.
          inset => Passed as an optional argument for situations where one wish to perturbate a data object besides this object's interal data.
                   I.E. you want to perturbate a parameter set that is not the parameters contained in this data object.
          setOnly => By default this function returns a data object of the same type as this object containing the new data set.  However,
                     there are times such as during playouts where one only wants the data set itself and not the data object. Setting this to true
                     returns the data set


        Out Signature
          newobj => Data Object of the same type as this Data Object that contains the newly perturbated Data Set.
           -- or if inset==True --
          newstruct => New Data Set created by this function without the Data Object wrapper
          
        '''
        OperationList = [self.mutate_sequence] 
        weights  = [1]

        #SatisfiesConstrains = False
        cnt = 0
        nodeDepth = node.getdepth()
        mutOperation = random.choices(OperationList,weights,k=1)[0]
        newPeptideDict = mutOperation(nodeDepth=nodeDepth)

        if setOnly:
            return newPeptideDict
        else:
            newobj = PeptideOptData(peptideDict=newPeptideDict)
            newobj.setevalutor(self.evaluator)
            return newobj
    #------------------------------------------------

    
    
    def newdataobject(self):
        '''
        Returns an initialized data object of the same type as this class. Used by the MCTS to create data objects
        for misc purposes and new node creations. 
        '''
        newobj = PeptideOptData(peptideDict=self.getPeptideDict())
        newobj.setevalutor(self.evaluator)
        return newobj

    
    
    
    def runsim(self, playouts=3, node=None, que_submit=que_submit):
        '''

        In Signature 
          playouts => Number of Random simulations to perform
        Out Signature
          structlist => Ordered list containing the actual data sets created by the random playouts.
          scorelist => Ordered list containing the matching objective scores 
        '''
        moves = 1
        peptideDictList=[] 
        energyList=[] 
        use_ML = False
        
        
        if node is not None: 
            print("\n\nNode Depth: %s"%(node.getdepth()))

        print('Playouts Value is:', playouts)

        if que_submit:
            energyList, peptideDictList = self._run_que_submit(playouts=playouts, node=node, n_ML_playouts=n_ML_playouts)
        else:
            energyList, peptideDictList = self._run_sequential(playouts=playouts, node=node, n_ML_playouts=n_ML_playouts)

        return energyList, peptideDictList
    #----------------------------------------------------



    def _run_sequential(self, playouts=3, node=None, n_ML_playouts=0):
        """
        Performs playouts Sequentially. n_ML_playouts is the number of ML-based playouts 
        """

        peptideDictList=[] 
        energyList=[]

        for playout in range(playouts):
            if playout < n_ML_playouts:
                # ML based selection
                moves = 5
                minscore = 0
                
                for iMove in range(moves):
                    try:
                        sampleObj = self.perturbate(node=node)
                    except: 
                        print('Failed during Perturbation with Data:', self.getPeptideDict())
                        raise

                    sampleDict = sampleObj.getPeptideDict()
                    pred_score, peptide = self.evaluator(sampleDict,
                                                label='{}00{}'.format(node.getid(),playout),
                                                use_ML=True)
                    if pred_score is None:
                        print('Not enough samples to train ML model')
                        newObj = sampleObj
                        break
                        
                    if pred_score < minscore:
                        newObj = sampleObj
                        minscore = pred_score
                print('ML based playout selection %s'%'-'.join(newObj.peptide))
            else:
                moves = 1
                # Random selection 
                for iMove in range(moves): ####
                    try:
                        newObj = self.perturbate(node=node)
                    except: 
                        print('Failed during Perturbation with Data:', self.getPeptideDict())
                        raise
                print('Random playout selection %s'%'-'.join(newObj.peptide))


            newDict = newObj.getPeptideDict()
                
            if node.getdepth() > 3: 
                reset=True
            else: 
                reset=False
                
            energy,finalPeptide = self.evaluator(newDict,
                                                label='{}00{}'.format(node.getid(),playout),
                                                reset=reset, use_ML=False)

            print("Playout %s Result: %s"%(playout, energy))
            energyList.append(energy)
            peptideDictList.append(finalPeptide)

        return energyList, peptideDictList
    #----------------------------------------------------


    def _run_que_submit(self, playouts=3, node=None, n_ML_playouts=0):
        """
        Performs playouts in parallel using job submission. Suitable for performing on clusters. 
        n_ML_playouts is the number of ML-based playouts 
        """
        joblist = []
        peptides = []
        for playout in range(playouts):
            if playout < n_ML_playouts:
                # ML based selection
                moves = 50
                minscore = 0
                
                for iMove in range(moves):
                    try:
                        sampleObj = self.perturbate(node=node)
                    except: 
                        print('Failed during Perturbation with Data:', self.getPeptideDict())
                        raise

                    sampleDict = sampleObj.getPeptideDict()
                    pred_score, peptide = self.evaluator(sampleDict,
                                                label='{}00{}'.format(node.getid(),playout),
                                                use_ML=True)

                    if pred_score < minscore:
                        newObj = sampleObj
                        minscore = pred_score
            else:
                moves = 1
                # Random selection 
                for iMove in range(moves): ####
                    try:
                        newObj = self.perturbate(node=node)
                    except: 
                        print('Failed during Perturbation with Data:', self.getPeptideDict())
                        raise


            newDict = newObj.getPeptideDict()
            peptides.append(newDict)
                
            if node.getdepth() > 3: 
                reset=True
            else: 
                reset=False

            # Submit job for the peptide
            job = que_playout_evaluate(newDict, queSystem)
            if job is not None:
                joblist.append(job)

        # Monitor all the jobs for completion
        queSystem.monitor(joblist)

        # Read Completed Peptides
        energyList, peptideDictList = que_playout_results(peptides)


        #    energy,finalPeptide = self.evaluator(newDict,
        #                                        label='{}00{}'.format(node.getid(),playout),
        #                                        reset=reset, use_ML=False)

        #    print("Playout %s Result: %s"%(playout, energy))
        #    energyList.append(energy)
        #    peptideDictList.append('{}00{}.{}'.format(node.getid(),playout,finalPeptide))

        return energyList, peptideDictList
#----------------------------------------------------            
    
    def computescore(self, que_submit=que_submit):
        '''
        Computes the objective function associated with the data contained within this object. Usually called when a new 
        node is created.
        '''
        peptideDict = self.getPeptideDict()
        #print(peptideDict)

        if que_submit:
            joblist = []

            # Submit job for the peptide
            job = que_playout_evaluate(peptideDict, queSystem)
            if job is not None:
                joblist.append(job)

            # Monitor all the jobs for completion
            queSystem.monitor(joblist)

            # Read Completed Peptides
            energyList, peptideDictList = que_playout_results([peptideDict])
            self.energy = energyList[0]

        else:
            self.energy,peptideDict = self.evaluator(peptideDict)
        return self.energy


    def getPeptidefp(self):
        # Fingerprint is not set. Read the PDB file get SMILES and then MorganFingerprints
        if self.fingerprint is None:     
            peptide = '-'.join(self.peptide)
            #m = Chem.MolFromPDBFile(pepStructuredir + '%s_aa.pdb'%peptide)
            m = Chem.rdmolfiles.MolFromSequence(self.short_sequence)

            if m is not None:
                self.smiles = Chem.MolToSmiles(m)
                fingerprint = GetMorganFingerprint(m, 3, useFeatures=True)
                return fingerprint
            else:
                print('Error Reading PDB file for %s'%peptide)
                return None
        return self.fingerprint

        
    def getuniqueness(self,node=None,nodelist=None,fingerprintfunc='mbtr'):
        '''
        Uniqueness Scoring Function. Returns a score based on how different this structure is compared to other
        sampled structures. This is designed to guide the exploration part of the MCTS selection formula to explore
        structures that are new. Failure to define this can result in bad optimization as the optimizer may tend toward
        finding degenerate Data Set combinations. 
        
        I.E. In structure optimization it's possible to have the same exact structure, but have
        it rotated or shift by a factor.  This function would contain a fingerprinting function that
        can distinguish between identical structures. 
        '''
        self.fingerprint = self.getPeptidefp()
        
        fingprntList = [node.data.getPeptidefp() for node in nodelist]
        distance = np.array([DataStructs.DiceSimilarity(self.fingerprint,fp) for fp in fingprntList])        
        uniqueness = 1 - np.median(distance)
        
        return uniqueness


    def setstructure(self, structPath):
        '''
          Sets the Data Set contained within this Data Object
        '''
        self.peptide = structPath.split('.')[-1].split('-')
        self.nblocks = len(self.peptide) 
        self.fingerprint = None
        self.smiles = None
        self.short_sequence = get_short_sequence(self.peptide)

    def datatostr(self):
        """
        Convert data object to string value to be writtin in file
        """
        return '-'.join(self.peptide)

    def convertstr(self, instr):
        """
        Convert string form of the data to be set as data object using setstructure module
        """
        return instr

    

