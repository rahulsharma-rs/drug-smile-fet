from setuptools import setup
import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name='dsfet',
    version='1.0',
    license='MIT',
    author='Rahul Sharma',
    author_email='rahul_rs.sharma@hotmail.com',
    long_description=long_description,
    description='This tool provides methods to extract meaningful features from drug SMILES for Machine Learning operation',
    url='https://github.com/rahulsharma-rs/drug-smile-fet',
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    keywords=['Drug SMILE', 'Feature Extraction', 'NLP'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
    #py_modules=['fe_1mol'],
    #package_dir={'dsfet'},
    install_requires = [
        'scikit-learn',
        'rdkit-pypi',
        'pandas'
    ]
)
