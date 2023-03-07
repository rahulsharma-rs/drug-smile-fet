import pandas as pd
from dsfet import fe_1mol
smile = 'N[C@](Br)(O)C'
lst = fe_1mol.getKmers(smile)
df = pd.read_csv('/Users/rahulsharma/Dropbox/UAB/drug-smile-fet/tests/SMILES_FeatureEngineered.csv')
df = df[['DRUG_NAME', 'PUBCHEM_ID', 'SMILES']].copy(deep=True)

df_1x = pd.read_csv('/Users/rahulsharma/Dropbox/UAB/drug-smile-fet/tests/SMILES_FeatureEngineered.csv')
df_1x.rename(columns={'Drug': 'DRUG_NAME', 'Cancer Type': 'TCGA_DESC'}, inplace=True)

#example of feature extraction from SMILES using NLP-based feature extraction method

#Train, Test, feature_sequences, feature_to_token_map = oneMolFeatureExtraction(trainSMILES=df, testSMILES=df_1x,ngram_list=[1,2,3,4,5,6,7,8])
Train, Test, feature_sequences, feature_to_token_map = fe_1mol.oneMolFeatureExtraction(trainSMILES=df, testSMILES=None,ngram_list=[1,2,3,4,5,6,7,8])

#example of feature extraction from SMILES using Morgan Fingerprints
result= fe_1mol.morganFingerprints(df, nBits=1024)
print('x')