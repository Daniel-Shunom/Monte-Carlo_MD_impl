# Import libraries
from copy import deepcopy
import pandas as pd

from joblib import Parallel, delayed
from tqdm import tqdm

from sklearn.model_selection import train_test_split
from sklearn.model_selection import GridSearchCV

import matplotlib.pyplot as plt
import sys

from sklearn.ensemble import RandomForestRegressor
import pickle

from os import path

from peptideUtilities import get_logP, get_norm_AP, get_short_sequence, get_smiles


class PeptideMLModel(object):
    '''
    ML model to make predictions for new peptides and possibly select cases with top scores
    '''
    #------------------------------------------------

    def __init__(self, model_type='RF',
                 datafname='/home/rbatra/Research/peptide_sequence/run_3mers/current_peptide_scores.csv',
                 model_fname = '/home/rbatra/Research/peptide_sequence/run_3mers/rf_AP_3mer.ml'):
        
        self.model_type = model_type
        self.datafname = datafname
        self.model = None
        self.Xcols = None
        self.ntrain = 0
        self.model_fname = model_fname
        self.pre_computed_fp = pd.read_csv('ml_model/fp_data.csv', index_col=0)

        if path.exists(model_fname):
            print('ML Model Found. Loading saved model')
            self.model = pickle.load(open(model_fname, 'rb'))
            self.Xcols = self.model.Xcols
            self.ntrain = self.model.ntrain

        
    def train(self, train_data_fname, fp_fname):
        
        rnseed = 1121
        ncvfold = 5
        scoring = 'neg_mean_absolute_error'
        njob = 24
        
        # Read Train Data File        
        if train_data_fname is not None:
            train_data = pd.read_csv(train_data_fname)
        else:
            train_data = pd.read_csv(self.datafname)

        
        # Get Fingerprint data
        # For old cases, fp taken from fp_fname file
        # For new cases, fp is computed and written in the fp_fname file

        #fp_data = pd.read_csv(fp_fname)

        # Start Model training
        fp_data = train_data.merge(self.pre_computed_fp, on=['Peptide','AP'], how='inner')
        print("Attempting ML training with %s points"%len(fp_data))

        Xcols = fp_data.columns[fp_data.columns.str.contains('fp_')]

        X = fp_data[Xcols]
        y = fp_data['AP']
        
        # Train RF model
        param_grid = {"n_estimators": range(100,101,100)}
        forest_reg = RandomForestRegressor(random_state=42, n_jobs=-1)

        # Perform Grid search for C and Gamma hyper-parameter        
        reg = GridSearchCV(forest_reg,
                           cv = ncvfold,
                           param_grid = param_grid,
                           scoring=scoring)

        reg.fit(X.values, y.values)

        self.model = deepcopy(reg.best_estimator_)
        self.Xcols = Xcols
        self.ntrain = len(fp_data)

        self.model.Xcols = Xcols
        self.model.ntrain = len(fp_data)
        
        pickle.dump(self.model, open(self.model_fname, 'wb'))
        
        print('ML model trained successfully with Ntrain %s\n\n' %self.ntrain)
        
        
    def predict_APscore(self, peptide):
        
        # Get Fingerprint data (To be updated later)
        #fp_data = pd.read_csv('ml_model/fp_data.csv', index_col=0)

        #fp_data = fp_data[['Peptide'] + list(fp_data.columns[fp_data.columns.str.contains('fp_')])]
        fp = self.pre_computed_fp.loc[self.pre_computed_fp['Peptide']=='-'.join(peptide),self.Xcols]

        '''
        # Fingerprint using fp code
        shortSeq = get_short_sequence(peptide)
        smiles = get_smiles(shortSeq)
        fp = pd.DataFrame([get_single_smile_pg_fp(smiles)])

        # Check is All Xcols are present
        newcols = self.Xcols[~self.Xcols.isin(fp.columns)]
        for newcol in newcols:
            fp[newcol] = 0

        fp = fp[self.Xcols].values
        '''

        # Make ML prediction
        pred_APscore = self.model.predict(fp)
        
        return pred_APscore
    
    def getScore(self, peptide):
        
        # Get hydrophobicity Score
        logP, norm_logP = get_logP('-'.join(peptide))
        
        # Get AP score
        try:
            AP = self.predict_APscore(peptide)
            norm_AP = get_norm_AP(AP)
        except:
            print("Failed with ML prediction")
            norm_AP = 0

            
        # Get overall score
        score = (norm_AP**2) * (norm_logP**0.5)
        #print('Predict Score: ', score)
        
        return score
