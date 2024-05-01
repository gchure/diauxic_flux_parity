import numpy as np 
import pkg_resources 
import pandas as pd
import frontmatter
import warnings
from io import StringIO
import pickle
import tqdm

__GENES__ = None
__STRAINS__ = None
def _load_colidict():
    global __GENES__
    f = pkg_resources.resource_filename('diaux', 'package_data/coli_gene_dict.pkl')
    with open(f, 'rb') as file:
        __GENES__ = pickle.load(file)

def _load_strains():
    global __STRAINS__
    f = pkg_resources.resource_filename('diaux','package_data/strain_database.pkl')
    with open(f, 'rb') as file:
        __STRAINS__ = pickle.load(file)

def standardize_strains(strains):
    """
    Standardizes strain names across different mutants given GC database 
    strain identifiers.

    Parameters
    ----------
    strains: str, list, or 1d numpy.ndarray
        List of GC Database strain identifiers (e.g. GC001, GC002) you wish to 
        standardize.

    Returns
    --------
    [shorthand, [genotype, lab_annotation], class] : list of lists
        A standardized list of lists with the standardized shorthand notations,
        the strain genotype, the annotation of the lab stock, and the strain
        class, which is either WT, Single KO, or Double KO.
    """
    global __STRAINS__ 
    if __STRAINS__ is None:
        _load_strains()
    
    # Peform the standardization
    out = [[], [[], []], []]
    for s in strains:
        sel = __STRAINS__[s]
        out[0].append(sel['shorthand'])
        out[1][0].append(sel['genotype'])
        out[1][1].append(sel['lab_annotation'])
        out[2].append(sel['class'])
    return out

def scrape_frontmatter(dirname, file='README.md'):
    """
    Reads the status of a given experimental dataset. This status is embedded
    in the README.md file as a YAML metadata block.

    Parameters
    ----------
    dirname : str
        Directory from which to parse.
    file: str
        Name of file containing YAML frontmatter. Default is 'README.md'

    Returns
    -------
    info : dict or pandas DataFrame
        A dictionary with all frontmatter keys and values.

    Raises
    ------
    UserWarning
        A UserWarning is raised if the scraped yaml frontmatter does not have
        a 'status' key or the value is not in `['accepted', 'rejected',
        'questionable']`.
    """
    # Grab file from directory.
    if dirname[-1] == '/':
        filename = '{}{}'.format(dirname, file)
    else:
        filename = '{}/{}'.format(dirname, file)

    # Scrape and return as desired.
    with open(filename) as f:
        info = frontmatter.load(f).metadata
        if 'status' in info.keys():
            info['status'] = info['status'].replace('\n', '').replace(' ', '')
        else:
            info = {}
    if 'status' not in info.keys():
        print(
            'Key `status` not found in metadata keys. Skipping {}'.format(dirname))
        info = {}
    elif info['status'] is None:
        print('Key `status` is missing. Skipping {}'.format(dirname))
        info = {}
    elif info['status'].lower() not in ['accepted', 'questionable', 'rejected']:
        print('Warning: Value `status: {}` not an acceptable flag. Skipping {}'.format(
            info['status'].lower(), dirname))
        info = {}
    return info