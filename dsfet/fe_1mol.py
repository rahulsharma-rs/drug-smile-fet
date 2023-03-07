import os, sys
from sklearn.preprocessing import StandardScaler, MaxAbsScaler, MinMaxScaler, QuantileTransformer, PowerTransformer,Normalizer
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.feature_extraction.text import CountVectorizer
from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np

from io import StringIO
sio = sys.stderr = StringIO()

two_char_mols = ['Cl', 'Br', 'He' , 'Li' , 'Be' , 'Ne' , 'Na' , 'Mg' ,
                 'Al' , 'Si' , 'Cl' , 'Ar' , 'Ca' , 'Sc' , 'Ti', 'Cr' ,
                 'Mn' , 'Fe' , 'Co' , 'Ni' , 'Cu' , 'Zn' , 'Ga' , 'Ge' ,
                 'As' , 'Se' , 'Br' , 'Kr' , 'Rb' , 'Sr' , 'Zr' , 'Nb' ,
                 'Mo' , 'Tc' , 'Ru' , 'Rh' , 'Pd' , 'Ag' , 'Cd' , 'In' ,
                 'Sn' , 'Sb' , 'Te' , 'Xe' , 'Cs' , 'Ba' , 'Hf' , 'Ta' ,
                 'Re' , 'Os' , 'Ir' , 'Pt' , 'Au' , 'Hg' , 'Tl' , 'Pb' ,
                 'Bi' , 'Po' , 'At' , 'Rn' , 'Fr' , 'Ra' , 'Rf' , 'Db' ,
                 'Sg' , 'Bh' , 'Hs' , 'Mt' , 'Ds' , 'Rg' , 'Cn' , 'Fl' ,
                 'Lv' , 'La' , 'Ce' , 'Pr' , 'Nd' , 'Pm' , 'Sm' , 'Eu' ,
                 'Gd' , 'Tb' , 'Dy' , 'Ho' , 'Er' , 'Tm' , 'Yb' , 'Lu' ,
                 'Ac' , 'Th' , 'Pa' , 'Np' , 'Pu' , 'Am' , 'Cm' , 'Bk' ,
                 'Cf' , 'Es' , 'Fm' , 'Md' , 'No' , 'Lr' ,'se','as']

# function to convert sequence strings into k-mer words, default size = 6 (hexamer words)
def getKmers1(sequence, size=6):
    if sequence != None:
        lst = []
        X = None
        for i in [1]:
            for x in range(len(sequence) - i + 1):
                try:

                    X = Chem.MolFromSmiles(sequence[x:x + i])
                    sio = sys.stderr = StringIO()
                except SyntaxError:
                    pass
                if X is None:
                    continue
                else:
                    lst.append(sequence[x:x + i])
        return lst
    else:
        return None

def getKmers(sequence, size=6):
    if sequence != None:
        lst = []
        X = None
        i = 2
        #for i in [2]:
        x = 0
        #for x in range(len(sequence) - i + 1):
        while x in range(len(sequence) - i + 1):
            if sequence[x:x + i] in two_char_mols:
                lst.append(sequence[x:x + i])
                x +=i
            else:
                for y in sequence[x:x + i]:
                    try:

                        X = Chem.MolFromSmiles(y)
                        sio = sys.stderr = StringIO()
                    except SyntaxError:
                        pass
                    if X is None:
                        continue
                    else:
                        lst.append(y)
                x += i
        return lst
    else:
        return None

def mapWords(lst=None,mdict=None):

    if lst!=None:
        tlst = []
        for x in lst:
            tlst.append(mdict[x])
        return tlst
    else:
        return None
"""
The function oneMolFeatureExtraction takes three arguments, trainSMILES, testSMILES, and ngram_list.

trainSMILES is a pandas DataFrame containing a column with SMILES strings that represent chemical structures. testSMILES is an optional argument that represents the same as trainSMILES, but for a separate set of test data.

ngram_list is a list of integers that specifies the range of n-grams to use for feature extraction.

The function performs the following operations:

It generates a list of all unique single character molecules from the SMILES strings in trainSMILES, and creates aliases for each molecule.
It creates a mapping between each single character molecule and its respective alias.
It maps the SMILES strings in trainSMILES to sequences of their corresponding aliases.
It performs feature extraction on the mapped sequences using a CountVectorizer, with n-gram sizes specified by the ngram_list argument.
It performs z-transformation on the extracted features using the StandardScaler.
It creates a new DataFrame containing the original trainSMILES DataFrame concatenated with the extracted features, and returns it.
If testSMILES is not None, it performs the same operations on the testSMILES DataFrame and returns a new DataFrame containing the concatenated extracted features.
In summary, this function takes a DataFrame with SMILES strings and performs feature extraction on them using a specified range of n-grams, returning a new DataFrame with the extracted features appended as additional columns. It also performs z-transformation on the extracted features for normalization. If test data is provided, it returns a separate DataFrame containing the extracted features for the test data. The function is useful for generating features for use in machine learning models that use chemical structures as input.
"""
def oneMolFeatureExtraction(trainSMILES=None,testSMILES=None,ngram_list =None):

    df = trainSMILES
    df_1x = testSMILES

    # Selecting single character-based legitimate molecule
    df['WORDS'] = df.apply(lambda x: getKmers(x['SMILES']), axis=1)
    if df_1x is not None:
        df_1x['WORDS'] = df_1x.apply(lambda x: getKmers(x['SMILES']), axis=1)

    mols = []  # all legit single character molecule for all drug SMILES
    for x in list(df['WORDS']):
        mols = mols + x

    setx = list(set(mols))  # Unique Single Character Molecule set
    molName = [f'mol{i}' for i in range(1, len(setx) + 1)]  # aliases for all single character Molecules
    molDict = {}  # Dictionary to map the Single Character Molecule to their respective aliases

    ctr = 1
    for i in setx:
        molDict.update({i: f'mol{ctr}'})
        ctr += 1
    mol = pd.DataFrame(np.asarray([setx, molName]).T, columns=['Molecule', 'tokenName'])
    #mol.to_csv('../../D2GNets/data/Molecule_token_map.csv', index=False)
    # words based on sequence of aliases of the drug SMILE
    df['WordMap']=df.apply(lambda x: mapWords(x['WORDS'],molDict), axis=1)
    if df_1x is not None:
        df_1x['WordMap']=df_1x.apply(lambda x: mapWords(x['WORDS'],molDict), axis=1)


    df_texts = list(df['WordMap']) # all aliases for all drug SMILES
    for item in range(len(df_texts)):
        df_texts[item] = ' '.join(df_texts[item])

    if df_1x is not None:
        df_texts1= list(df_1x['WordMap']) # all aliases for all drug SMILES
        for item in range(len(df_texts1)):
            df_texts1[item] = ' '.join(df_texts1[item])


    n_gram_list = ngram_list
    ctr =0
    # Creating the Bag of Words model using CountVectorizer()
    # This is equivalent to k-mer counting
    # The n-gram size of 4 was previously determined by testing
    for ind in n_gram_list:  #looping to create 8 type of n-gram vectors
        #feature extraction based on Counts
        cv = CountVectorizer(ngram_range=(ind,ind))
        X = cv.fit_transform(df_texts)
        if df_1x is not None:
            X1 = cv.transform(df_texts1)
        #Conveting the extracted features to thier term frequencies
        tf_transformer = TfidfTransformer(use_idf=False).fit(X)
        X = tf_transformer.transform(X)
        if df_1x is not None:
            X1 = tf_transformer.transform(X1)
        count_vect_df = pd.DataFrame(X.todense(), columns=cv.vocabulary_.keys())
        if df_1x is not None:
            count_vect_df1 = pd.DataFrame(X1.todense(), columns=cv.vocabulary_.keys())
        if ctr ==0:
            ctr+=1
            dff = count_vect_df
            #dff = pd.concat([df.reset_index(), count_vect_df], axis =1, ignore_index= False)
            if df_1x is not None:
                dff_1x = count_vect_df1
                #dff_1x = pd.concat([df_1x.reset_index(), count_vect_df1], axis =1, ignore_index= False)
        else:
            dff = pd.concat([dff, count_vect_df], axis =1, ignore_index= False)
            if df_1x is not None:
                dff_1x = pd.concat([dff_1x, count_vect_df1], axis =1, ignore_index= False)

    print('number of drug selected: %d'%dff.shape[0])
    print('number of features created: %d'%dff.shape[1])

    #saving the dataFrame for future reference
    #dff.to_csv('../../D2GNets/data/SMILES_FeatureEngineered_new.csv',index=False)

    #dff_1x.to_csv('../../D2GNets/data/SMILES_FeatureEngineered_Repurposing_drugs.csv',index=False)

    #Z-transforamtion of each feature of NLP-extracted drug features
    max_abs_scaler = StandardScaler()
    d1_1 = max_abs_scaler.fit_transform(dff[dff.columns])
    d1 = pd.DataFrame(d1_1,columns=dff.columns.to_list())
    #dff1 = pd.concat([dff[dff.columns],d1],axis=1)
    #dff1= dff.copy(deep=True)
    if df_1x is not None:
        d1_1x = max_abs_scaler.fit_transform(dff_1x[dff_1x.columns])
        d1_1x = pd.DataFrame(d1_1,columns=dff_1x.columns.to_list())
        #dff1_1x = pd.concat([dff_1x[dff_1x.columns],d1],axis=1)

    extracted_features = []
    for x in dff.columns.to_list():
        t_list = x.split(sep=' ')
        temp = []
        for y in t_list:
            temp.append((mol[mol['tokenName']==y])['Molecule'].values[0])
        extracted_features.append(' '.join(temp))

    trainSMILEdf = pd.concat([df.reset_index(), d1], axis =1, ignore_index= False)
    if df_1x is not None:
        testSMILESdf = pd.concat([df_1x.reset_index(), d1_1x], axis =1, ignore_index= False)
    if df_1x is not None:
        return trainSMILEdf, testSMILESdf, extracted_features, mol
    if df_1x is None:
        return trainSMILEdf, None, extracted_features, mol


"""
The function morganFingerPrint takes two arguments, data and nBits.

data is a pandas DataFrame containing a column with SMILES strings that represent chemical structures.

nBits is an integer that specifies the size of the fingerprint vector to be generated.

The function iterates over each row of the data DataFrame, converts the SMILES string to a molecular object using the RDKit library, generates the Morgan fingerprint for the molecule using the GetMorganFingerprintAsBitVect function, and appends the fingerprint as a list to a list of fingerprints called lst.

The list of fingerprints is then converted to a pandas DataFrame df, with each fingerprint being a row in the DataFrame. The column names of the DataFrame are labeled F1, F2, ..., FnBits, where nBits is the number of bits in the fingerprint vector.

Finally, the original data DataFrame is concatenated with the new df DataFrame, and the resulting DataFrame is returned as the output of the function.

In summary, this function takes a DataFrame with SMILES strings and generates Morgan fingerprints for each molecule, returning a new DataFrame with the fingerprints appended as additional columns.
"""
def morganFingerPrint(data=None, nBits=512):
    #smile = 'N[C@](Br)(O)C'
    nbts=nBits
    lst = []
    df = pd.DataFrame()
    if data is None:
        print('No data provided')
        return None
    else:
        for index, row in data.iterrows():
            ms = Chem.MolFromSmiles(row['SMILES'])
            pairFps = AllChem.GetMorganFingerprintAsBitVect(ms,2,nBits=nbts)
            lst.append(pairFps.ToList())
        df = pd.DataFrame(lst)
        df.columns = ['F'+str(x) for x in range(1,nbts+1)]
        df = pd.concat([data.reset_index(), df], axis =1, ignore_index= False)
    return df
'''
smile = 'N[C@](Br)(O)C'
lst = getKmers1(smile)
df = pd.read_csv('/Users/rahulsharma/Dropbox/UAB/drug-smile-fet/tests/Drugs_for_repurposing.csv')
df = df[['DRUG_NAME', 'PUBCHEM_ID', 'SMILES']].copy(deep=True)

#df_1x = pd.read_csv('/Users/rahulsharma/Dropbox/UAB/drug-smile-fet/tests/Drugs_for_repurposing.csv')
#df_1x.rename(columns={'Drug': 'DRUG_NAME', 'Cancer Type': 'TCGA_DESC'}, inplace=True)
#Train, Test, feature_sequences, feature_to_token_map = oneMolFeatureExtraction(trainSMILES=df, testSMILES=df_1x,ngram_list=[1,2,3,4,5,6,7,8])
#Train, Test, feature_sequences, feature_to_token_map = oneMolFeatureExtraction(trainSMILES=df, testSMILES=None,ngram_list=[1,2,3,4,5,6,7,8])
train = morganFingerPrint(data=df, nBits=512)
print('x')
'''

