import os 
import pandas as pd
import subprocess
from peptideMLObject import PeptideMLModel
from peptideUtilities import get_logP, get_norm_AP
from time import sleep

workdir='./'

ml_train_freq = 100
model = PeptideMLModel(datafname='./data/current_peptide_scores.csv',
                        model_fname = './ml_model/rf_AP_3mer.ml')


def runGROMACS(peptide, workdir=workdir):
    
    mcts_path = os.getcwd()
    print(workdir, mcts_path)
    os.chdir(workdir)

    # Run GROMACS script
    command = ['./run_all.sh'] + peptide.split('-')
    out = subprocess.run(command)

    try:
        # Collect Results
        with open('5_Analysis/%s_sasa_init.xvg'%peptide, 'r') as f:
            last_line = f.readlines()[-1]
            f.close()
        init_sasa = float(last_line.split()[-1])

        with open('5_Analysis/%s_sasa_end.xvg'%peptide, 'r') as f:
            last_line = f.readlines()[-1]
            f.close()
        end_sasa = float(last_line.split()[-1])
        
        AP = init_sasa/end_sasa
        norm_AP = get_norm_AP(AP)
        logP, norm_logP = get_logP(peptide)

        score = (norm_AP**2) * (norm_logP**0.5)

    except:
        print("GROMACS Failed for peptide: ", peptide)  
        AP, logP, norm_AP, norm_logP, score = (-1,-1,-1,-1,-1)


    # Return to original Path
    os.chdir(mcts_path)
    
    return [AP, logP, norm_AP, norm_logP, score], peptide




# Evaluator for Peptide Work
def peptideEvaluate(newDict,label='0',minimize=False,reset=False, use_ML=False):

    for trial in range(5):
        try:
            peptide_info = pd.read_csv(workdir + 'data/peptide_scores.csv')
            this_run = pd.read_csv('data/current_peptide_scores.csv')
            break
        except:
            print('Error Reading Peptide Info File. Trying Again')
            sleep(2)
            continue

    if use_ML:
        # Check for first time training
        len_train_data = len(this_run)

        #if (len_train_data > 200) and (model.model is None):
        #    model.train(train_data_fname=workdir + 'data/current_peptide_scores.csv', fp_fname=workdir + 'ml_model/fp_data.csv')

        # Check if ML model (re)training is required
        if (len_train_data - model.ntrain) > ml_train_freq:
            print("\n\n(Re)Training ML model")
            model.train(train_data_fname=workdir + 'data/current_peptide_scores.csv', fp_fname=workdir + 'ml_model/fp_data.csv')

        if model.model is None:
            return None, None


        pred_score = model.getScore(newDict['peptide'])

        return -pred_score, '-'.join(newDict['peptide'])  # Negative of score since we want to maximize


    peptide = '-'.join(newDict['peptide'])
    #peptide_info = pd.read_csv('/Users/rbatra/Research/Research_ANL/Peptide_assembly/pastWork/triPepData_v2.csv')

    if peptide in peptide_info['Peptide'].values:
        print('%s MD ran previously, fetching MD evaluated score from file'%peptide)
        score = peptide_info.loc[peptide_info['Peptide']==peptide, 'Score'].values
        assert len(score)==1


        # Add information of current run
        if peptide not in this_run['Peptide'].values:
            this_run = this_run.append(peptide_info.loc[peptide_info['Peptide'] == peptide])
            this_run.to_csv('data/current_peptide_scores.csv', index=False)

        #new_info = peptide_info.loc[peptide_info['Peptide']==peptide]
        #assert len(new_info)==1

        #peptide_info_current = pd.read_csv(workdir + 'peptide_scores.csv')
        #if len(peptide_info_current.loc[peptide_info_current['Peptide']==peptide]) == 0:
        #    peptide_info_current = peptide_info_current.append(new_info, ignore_index=True)
        #    peptide_info_current.to_csv(workdir + 'peptide_scores.csv', index=False)

        return -score[0], peptide  # Negative of score since we want to maximize

    else:
        print('Running GROMACS for: ', peptide)

        scores, peptide = runGROMACS(peptide, workdir=workdir)
        #addPeptideInfo(scores, peptide)
        
        if scores[-1] != -1:            
            new_info = pd.DataFrame([{'Peptide': peptide,
                 'AP': scores[0],
                 'logP': scores[1],
                 'norm_AP': scores[2],
                 'norm_logP': scores[3],                 
                 'Score': scores[4]}])

            peptide_info = peptide_info.append(new_info, ignore_index=True)
            peptide_info.to_csv(workdir + 'data/current_peptide_scores.csv', index=False)
            
        return -scores[-1], peptide


def que_playout_evaluate(newDict, queSystem):

    peptide_info = pd.read_csv(workdir + 'data/current_peptide_scores.csv')
    peptide = '-'.join(newDict['peptide'])

    if peptide in peptide_info['Peptide'].values:
        print('Ran previously, Fetching Score for: ', peptide)
        score = peptide_info.loc[peptide_info['Peptide']==peptide, 'Score'].values
        assert len(score)==1

        return None

    mcts_path = os.getcwd()
    os.chdir(workdir)

    with open("submit.sh", "rt") as fin:
        with open("tmp_submit.sh", "wt") as fout:
            for line in fin:
                fout.write(line.replace('PepXXX', ' '.join(newDict['peptide'])))

    job = queSystem.createjob("tmp_submit.sh")
    os.chdir(mcts_path)

    return job

def que_playout_results(peptides):

    energyList = []
    peptideList = []

    mcts_path = os.getcwd()
    os.chdir(workdir)

    peptide_info = pd.read_csv(workdir + 'data/current_peptide_scores.csv')
    for newDict in peptides:
        peptide = '-'.join(newDict['peptide'])

        if peptide in peptide_info['Peptide'].values:
            print('Ran previously, Fetching Score for: ', peptide)
            score = peptide_info.loc[peptide_info['Peptide']==peptide, 'Score'].values
            assert len(score)==1

            energyList.append(-score[0])
            peptideList.append(peptide)

        else:

            try:
                # Collect Results
                with open('5_Analysis/%s_sasa_init.xvg'%peptide, 'r') as f:
                    last_line = f.readlines()[-1]
                    f.close()
                init_sasa = float(last_line.split()[-1])

                with open('5_Analysis/%s_sasa_end.xvg'%peptide, 'r') as f:
                    last_line = f.readlines()[-1]
                    f.close()
                end_sasa = float(last_line.split()[-1])
                
                AP = init_sasa/end_sasa
                norm_AP = get_norm_AP(AP)
                logP, norm_logP = get_logP(peptide)

                score = (norm_AP**2) * (norm_logP**0.5) 

            except:
                print("GROMACS Failed for peptide: ", peptide)  
                AP, logP, norm_AP, norm_logP, score = (-1,-1,-1,-1,-1)            

            # Add score to the Overall csv file
            if score != -1:            
                new_info = pd.DataFrame([{'Peptide': peptide,
                     'AP': AP,
                     'logP': logP,
                     'norm_AP': norm_AP,
                     'norm_logP': norm_logP,
                     'Score': score}])

                peptide_info = peptide_info.append(new_info, ignore_index=True)
            
            energyList.append(-score)
            peptideList.append(peptide)

    peptide_info.to_csv(workdir + 'data/current_peptide_scores.csv', index=False)
    os.chdir(mcts_path)

    return energyList, peptideList

