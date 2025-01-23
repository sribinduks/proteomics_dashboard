# File: data_loader.py

import os
import pandas as pd
from glob import glob
import math
import requests
from dotenv import load_dotenv
load_dotenv()

X_CDD_token=os.environ.get('X_CDD_TOKEN')
image_url = 'http://158.101.124.232:3060/image?smiles='

def get_smile_url(words):
    if words:
        vault_id = 5106
        base_url = f"https://app.collaborativedrug.com/api/v1/vaults/{vault_id}/molecules?names="
        headers = {'X-CDD-token':X_CDD_token}
        url = base_url + words
        response = requests.request("Get",url,headers=headers)
        if response.status_code == 200:
            response = response.json()
            objects = response['objects']
            smiles = objects[0]['smiles'].replace('#','%23').replace('+','%2B')
            return image_url + smiles
        else:
            return None
    else:
        return None
    
    
def load_data(folder):
    files = glob(os.path.join(folder,'ip*_output_CPD*.csv'))
    data = {}
    for file in files:
        # print(file)
        df = pd.read_csv(file)
        if 'P.Value' in df.columns.tolist():
            df['-logp']=df['P.Value'].apply(lambda x: -math.log10(x))
        file = os.path.basename(file)
        data[file] = df
    return data


def load_txt(folder):
    files = glob(os.path.join(folder,'*_PeptideGroups.txt'))
    data = {}
    for file in files:
        df = pd.read_csv(file, delimiter='\t')
        file = os.path.basename(file)
        data[file] = df
    return data

def get_compounds(data):
    compounds = dict()
    for file_name in data:
        if "CPD." in file_name:
            compound = file_name.split('output_')[1].split('_vs')[0]
            cpmd_id = compound.split('.')[1].split('_')[0]
            # add 0 before cmpd_id to reach 7 digits
            cpmd_id = cpmd_id.zfill(7)
            compound = 'CPD-' + cpmd_id
            compounds[file_name] = compound
        else:
            pass
    return compounds

def get_experiment(txt_data):
    experiments = {}
    for file_name in txt_data:
        experiment = file_name.rsplit('_',2)[0]
        experiments[file_name] = experiment
    return experiments

def get_txtfile_from_data(datafile,txt_data):
    fn_pre = '_'.join(datafile.split('_')[:3])
    for file_name in txt_data:
        if fn_pre in file_name:
            return [file_name]
