from setuptools import setup
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='drug-smile-fet',
    version='0.2',
    packages=['dsfet'],
    url='https://github.com/rahulsharma-rs/drug-smile-fet',
    license='MIT',
    author='Rahul Sharma',
    author_email='rahul_rs.sharma@hotmail.com',
    long_description=long_description,
    description='This tool provides methods to extract meaningful features from drug SMILES for Machine Learning operation',
    long_description_content_type="text/markdown",
    keywords=['Drug SMILE', 'Feature Extraction', 'NLP'],
    #description='This tool provides methods to extract meaningful features from drug SMILES for Machine Learning operation'
)
