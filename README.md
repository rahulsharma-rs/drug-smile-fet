## Pre-requisites
Install RdKit library:
- instllation through [anaconda](https://anaconda.org/rdkit/rdkit)
- instllation through [pip](https://pypi.org/project/rdkit-pypi/) 

## Usage

- Make sure you have Python installed in your system.
- Run Following command in the CMD.
 ```
  pip install drug-smile-fet
  ```
## Example

 ```
# example.py
from dsfet import fe_1mol
import pandas as pd
train_smiles = {'DRUG_NAME': {0: 'Luminespib', 1: 'Trametinib', 2: 'Venetoclax', 3: 'Olaparib', 4: 'Axitinib'},
                'PUBCHEM_ID': {0: 135539077.0, 1: 11707110.0, 2: 49846579.0, 3: 23725625.0, 4: 6450551.0},
                'SMILES': {0: 'CCNC(=O)C1=NOC(=C1C2=CC=C(C=C2)CN3CCOCC3)C4=CC(=C(C=C4O)O)C(C)C',
                           1: 'CC1=C2C(=C(N(C1=O)C)NC3=C(C=C(C=C3)I)F)C(=O)N(C(=O)N2C4=CC=CC(=C4)NC(=O)C)C5CC5',
                           2: 'CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)C',
                           3: 'C1CC1C(=O)N2CCN(CC2)C(=O)C3=C(C=CC(=C3)CC4=NNC(=O)C5=CC=CC=C54)F',
                           4: 'CNC(=O)C1=CC=CC=C1SC2=CC3=C(C=C2)C(=NN3)/C=C/C4=CC=CC=N4'}
                }
train_smiles_df = pd.DataFrame(data=train_smiles)

test_smile = train_smiles
test_smile_df = pd.DataFrame(test_smile)

#Example 1: to call NLP-based feature extraction method
#Train, Test, feature_sequences, feature_to_token_map = fe_1mol.oneMolFeatureExtraction(trainSMILES=train_smiles_df, testSMILES=train_smiles_df,ngram_list=[1,2,3,4,5,6,7,8])
Train, Test, feature_sequences, feature_to_token_map = fe_1mol.oneMolFeatureExtraction(trainSMILES=train_smiles_df, testSMILES=None,ngram_list=[1,2,3,4,5,6,7,8])

#Example 2: to call Morgan Fingerprints based feature extraction method
#nBits is the number of bits in the fingerprint
result= fe_1mol.morganFingerPrint(train_smiles_df, nBits=1024)

```
### Note: 
The input to the method ```oneMolFeatureExtraction() and morganFingerPrint()``` must be a pandas DataFrame and the Drug SMILES column name must be in uppercase:
- e.g., SMILES

### Cite us at:
Rahul Sharma, Ehsan Saghapour, Jake Y. Chen, An NLP-based Technique to Extract Meaningful Features from Drug SMILES, iScience, 2024, 109127, ISSN 2589-0042, https://doi.org/10.1016/j.isci.2024.109127.
(https://www.sciencedirect.com/science/article/pii/S2589004224003481)
### Bibtex reference
@article{sharma2024nlp,
 
  title={An NLP-based Technique to Extract Meaningful Features from Drug SMILES},
  
  author={Sharma, Rahul and Saghapour, Ehsan and Chen, Jake Y},
  
  journal={iScience},
  
  year={2024},
  
  publisher={Elsevier}

}


