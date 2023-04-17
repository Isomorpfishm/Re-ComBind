from __future__ import print_function
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import oddt
from oddt import interactions

from rdkit import Chem, DataStructs
from rdkit.Chem import MACCSkeys
from rdkit.Chem.rdmolfiles import MolFromPDBFile
from rdkit.Chem.rdmolfiles import MolFromSmarts
from rdkit.Chem.rdmolops import AddHs
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdFMCS import FindMCS
from rdkit.Chem.rdMolAlign import CalcRMS

from spyrmsd import io, rmsd

import pandas as pd
import numpy as np
from numpy.polynomial import Polynomial

import scipy
from scipy.optimize import curve_fit
from scipy.stats.stats import pearsonr
from scipy.spatial.distance import cdist
from scipy.stats import gaussian_kde as kde
from sklearn.preprocessing import MinMaxScaler


import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib

from itertools import product
import os
import glob
import math



scaler = MinMaxScaler()

settings = {}
settings['nonpolar'] = {6:1.7, 9:1.47, 17:1.75, 35:1.85, 53:1.98}
settings['hbond_dist_cut'] = 3.5
settings['hbond_angle_cut'] = 180.0
settings['hbond_dist_bound'] = 2.5
settings['hbond_angle_bound'] = 90.0

settings['hphob_dist_cut'] = 4.0
settings['hphob_dist_bound'] = 1.25

settings['contact_scale_cut'] = 1.75
settings['contact_scale_opt'] = 1.25

settings['sbridge_dist_cut'] = 5.0
settings['sbridge_dist_bound'] = 3.25
settings['saltbridge_resonance'] = True

settings['pistack_dist_cut'] = 6.5
settings['pistack_dist_bound'] = 3.8
settings['pistack_angle_bound'] = 60.0
settings['pistack_angle_cut'] = 90.0

settings['pication_dist_cut'] = 6.5
settings['pication_dist_bound'] = 4.3
settings['pication_angle_cut'] = 90.0
settings['pication_angle_bound'] = 30.0

settings['xbond_dist_cut'] = 5.0
settings['xbond_angle_cut'] = 150.0
settings['xbond_dist_bound'] = 3.5
settings['xbond_angle_bound'] = 135.0

BASE_ANGLES = np.array((0, 180, 120, 109.5, 90), dtype=float)

def angle(p1, p2, p3):
    """Returns an angle from a series of 3 points (point #2 is centroid).
    Angle is returned in degrees.
    Parameters
    ----------
    p1,p2,p3 : numpy arrays, shape = [n_points, n_dimensions]
        Triplets of points in n-dimensional space, aligned in rows.
    Returns
    -------
    angles : numpy array, shape = [n_points]
        Series of angles in degrees
    """
    v1 = p1 - p2
    v2 = p3 - p2
    return angle_2v(v1, v2)


def angle_2v(v1, v2):
    """Returns an angle between two vecors.Angle is returned in degrees.
    Parameters
    ----------
    v1,v2 : numpy arrays, shape = [n_vectors, n_dimensions]
        Pairs of vectors in n-dimensional space, aligned in rows.
    Returns
    -------
    angles : numpy array, shape = [n_vectors]
        Series of angles in degrees
    """
    # better than np.dot(v1, v2), multiple vectors can be applied
    dot = (v1 * v2).sum(axis=-1)
    norm = np.linalg.norm(v1, axis=-1) * np.linalg.norm(v2, axis=-1)
    return np.degrees(np.arccos(np.clip(dot/norm, -1, 1)))
