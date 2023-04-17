import numpy as np
import pandas as pd

import oddt
from oddt import interactions
from oddt.spatial import distance

from utils import settings
from utils import angle, angle_2v

import math
import glob



#######################################################################################
################################### Initialisation ####################################
#######################################################################################

def hbonds_oddt(protein, ligand, _int=True, cutoff=settings['hbond_dist_cut'], tolerance=30, mol1_exact=False, mol2_exact=False):
    a1d2 = interactions.hbond_acceptor_donor(protein, ligand, cutoff, tolerance, mol2_exact)
    a2d1 = interactions.hbond_acceptor_donor(ligand, protein, cutoff, tolerance, mol1_exact)
    
    dist_a1d2, dist_a2d1 = [], []
    a_resname_a1d2, d_resname_a1d2 = [], []
    a_resname_a2d1, d_resname_a2d1 = [], []
    a_resnum_a1d2, d_resnum_a1d2 = [], []
    a_resnum_a2d1, d_resnum_a2d1 = [], []
    a_atom_a1d2, d_atom_a1d2 = [], []
    a_atom_a2d1, d_atom_a2d1 = [], []
    angle_a1d2, angle_a2d1 = [], []
    type_a1d2, type_a2d1 = [], []
    
    """
    a1d2 : protein as acceptor, ligand as donor
    """
    
    angle1 = angle(a1d2[1]['coords'][:, np.newaxis, :],
                   a1d2[0]['coords'][:, np.newaxis, :],
                   a1d2[0]['neighbors'])
    angle2 = angle(a1d2[0]['coords'][:, np.newaxis, :],
                   a1d2[1]['coords'][:, np.newaxis, :],
                   a1d2[1]['neighbors'])
    
    for i in range(len(a1d2[0])):
        dist = math.dist(a1d2[0]['coords'][i], a1d2[1]['coords'][i])
        dist_a1d2.append(dist)
        a_resname, d_resname = a1d2[0]['resname'][i], a1d2[1]['resname'][i]
        a_resnum, d_resnum = a1d2[0]['resnum'][i], a1d2[1]['resnum'][i]
        a_atom, d_atom = a1d2[0]['atomtype'][i], a1d2[1]['atomtype'][i]
        
        if _int == True:
            if a_resname == '':
                a_resname = 'LIG'
            elif d_resname == '':
              d_resname = 'LIG'
        else:
            pass
            
        a_resname_a1d2.append(a_resname)
        d_resname_a1d2.append(d_resname)
        
        a_resnum_a1d2.append(a_resnum)
        d_resnum_a1d2.append(d_resnum)
        
        a_atom_a1d2.append(a_atom)
        d_atom_a1d2.append(d_atom)
        
        angle_a1d2.append(angle2[i][0])
        type_a1d2.append("hbond_donor")
        
    """    
    a2d1 : protein as donor, ligand as acceptor
    """ 
        
    angle1 = angle(a2d1[1]['coords'][:, np.newaxis, :],
                   a2d1[0]['coords'][:, np.newaxis, :],
                   a2d1[0]['neighbors'])
    angle2 = angle(a2d1[0]['coords'][:, np.newaxis, :],
                   a2d1[1]['coords'][:, np.newaxis, :],
                   a2d1[1]['neighbors'])
        
    for i in range(len(a2d1[0])):
        dist = math.dist(a2d1[0]['coords'][i], a2d1[1]['coords'][i])
        dist_a2d1.append(dist)
        a_resname, d_resname = a2d1[0]['resname'][i], a2d1[1]['resname'][i]
        a_resnum, d_resnum = a2d1[0]['resnum'][i], a2d1[1]['resnum'][i]
        a_atom, d_atom = a2d1[0]['atomtype'][i], a2d1[1]['atomtype'][i]
        
        if _int == True:
            if a_resname == '':
                a_resname = 'LIG'
            elif d_resname == '':
                d_resname = 'LIG'
        else:
            pass
            
        a_resname_a2d1.append(a_resname)
        d_resname_a2d1.append(d_resname)
        
        a_resnum_a2d1.append(a_resnum)
        d_resnum_a2d1.append(d_resnum)

        a_atom_a2d1.append(a_atom)
        d_atom_a2d1.append(d_atom)
                
        angle_a2d1.append(angle2[i][0])
        type_a2d1.append("hbond_acceptor")
        
    df = pd.DataFrame({"type"             : type_a1d2      + type_a2d1     ,
                       "acceptor_resname" : a_resname_a1d2 + a_resname_a2d1,
                       "acceptor_resnum"  : a_resnum_a1d2  + a_resnum_a2d1 ,
                       "acceptor_atom"    : a_atom_a1d2    + a_atom_a2d1   ,
                       "donor_resname"    : d_resname_a1d2 + d_resname_a2d1,
                       "donor_resnum"     : d_resnum_a1d2  + d_resnum_a2d1 ,
                       "donor_atom"       : d_atom_a1d2    + d_atom_a2d1   ,
                       "distance"         : dist_a1d2      + dist_a2d1     ,
                       "angle"            : angle_a1d2     + angle_a2d1    ,
                       "ring1_resname"    : [0 for i in range(len(type_a1d2 + type_a2d1))],
                       "ring1_resnum"     : [0 for i in range(len(type_a1d2 + type_a2d1))],
                       "ring2_resname"    : [0 for i in range(len(type_a1d2 + type_a2d1))],
                       "ring2_resnum"     : [0 for i in range(len(type_a1d2 + type_a2d1))],
                       "cation_atom"      : [0 for i in range(len(type_a1d2 + type_a2d1))],
                       })

    return df


def xbonds_oddt(protein, ligand, _int=True, cutoff=settings['xbond_dist_cut'], tolerance=30, mol1_exact=False, mol2_exact=False):
    a1h2 = interactions.halogenbond_acceptor_halogen(protein, ligand, cutoff, tolerance, mol2_exact)
    a2h1 = interactions.halogenbond_acceptor_halogen(ligand, protein, cutoff, tolerance, mol1_exact)
    
    dist_a1h2, dist_a2h1 = [], []
    a_resname_a1h2, h_resname_a1h2 = [], []
    a_resname_a2h1, h_resname_a2h1 = [], []
    a_resnum_a1h2, h_resnum_a1h2 = [], []
    a_resnum_a2h1, h_resnum_a2h1 = [], []
    a_atom_a1h2, h_atom_a1h2 = [], []
    a_atom_a2h1, h_atom_a2h1 = [], []
    angle_a1h2, angle_a2h1 = [], []
    type_a1h2, type_a2h1 = [], []
    
    """
    a1h2 : protein as acceptor, ligand as halogen donor
    """
    
    angle1 = angle(a1h2[1]['coords'][:, np.newaxis, :],
                   a1h2[0]['coords'][:, np.newaxis, :],
                   a1h2[0]['neighbors'])
    angle2 = angle(a1h2[0]['coords'][:, np.newaxis, :],
                   a1h2[1]['coords'][:, np.newaxis, :],
                   a1h2[1]['neighbors'])
    
    for i in range(len(a1h2[0])):
        dist = math.dist(a1h2[0]['coords'][i], a1h2[1]['coords'][i])
        dist_a1h2.append(dist)
        a_resname, h_resname = a1h2[0]['resname'][i], a1h2[1]['resname'][i]
        a_resnum, h_resnum = a1h2[0]['resnum'][i], a1h2[1]['resnum'][i]
        a_atom, h_atom = a1h2[0]['atomtype'][i], a1h2[1]['atomtype'][i]
        
        if _int == True:
            if a_resname == '':
                a_resname = 'LIG'
            elif h_resname == '':
              h_resname = 'LIG'
        else:
            pass
            
        a_resname_a1h2.append(a_resname)
        h_resname_a1h2.append(h_resname)
        
        a_resnum_a1h2.append(a_resnum)
        h_resnum_a1h2.append(h_resnum)
        
        a_atom_a1h2.append(a_atom)
        h_atom_a1h2.append(h_atom)
        
        angle_a1h2.append(angle2[i][0])
        type_a1h2.append("xbond_donor")
        
    """    
    a2h1 : protein as donor, ligand as halogen acceptor
    """ 
        
    angle1 = angle(a2h1[1]['coords'][:, np.newaxis, :],
                   a2h1[0]['coords'][:, np.newaxis, :],
                   a2h1[0]['neighbors'])
    angle2 = angle(a2h1[0]['coords'][:, np.newaxis, :],
                   a2h1[1]['coords'][:, np.newaxis, :],
                   a2h1[1]['neighbors'])
        
    for i in range(len(a2h1[0])):
        dist = math.dist(a2h1[0]['coords'][i], a2h1[1]['coords'][i])
        dist_a2h1.append(dist)
        a_resname, h_resname = a2h1[0]['resname'][i], a2h1[1]['resname'][i]
        a_resnum, h_resnum = a2h1[0]['resnum'][i], a2h1[1]['resnum'][i]
        a_atom, h_atom = a2h1[0]['atomtype'][i], a2h1[1]['atomtype'][i]
        
        if _int == True:
            if a_resname == '':
                a_resname = 'LIG'
            elif h_resname == '':
                h_resname = 'LIG'
        else:
            pass
            
        a_resname_a2h1.append(a_resname)
        h_resname_a2h1.append(h_resname)
        
        a_resnum_a2h1.append(a_resnum)
        h_resnum_a2h1.append(h_resnum)

        a_atom_a2h1.append(a_atom)
        h_atom_a2h1.append(h_atom)
                
        angle_a2h1.append(angle2[i][0])
        type_a2h1.append("xbond_acceptor")
        
    df = pd.DataFrame({"type"             : type_a1h2      + type_a2h1     ,
                       "acceptor_resname" : a_resname_a1h2 + a_resname_a2h1,
                       "acceptor_resnum"  : a_resnum_a1h2  + a_resnum_a2h1 ,
                       "acceptor_atom"    : a_atom_a1h2    + a_atom_a2h1   ,
                       "donor_resname"    : h_resname_a1h2 + h_resname_a2h1,
                       "donor_resnum"     : h_resnum_a1h2  + h_resnum_a2h1 ,
                       "donor_atom"       : h_atom_a1h2    + h_atom_a2h1   ,
                       "distance"         : dist_a1h2      + dist_a2h1     ,
                       "angle"            : angle_a1h2     + angle_a2h1    ,
                       "ring1_resname"    : [0 for i in range(len(type_a1h2 + type_a2h1))],
                       "ring1_resnum"     : [0 for i in range(len(type_a1h2 + type_a2h1))],
                       "ring2_resname"    : [0 for i in range(len(type_a1h2 + type_a2h1))],
                       "ring2_resnum"     : [0 for i in range(len(type_a1h2 + type_a2h1))],
                       "cation_atom"      : [0 for i in range(len(type_a1h2 + type_a2h1))],
                       })

    return df


def hphob_oddt(protein, ligand, _int=True, cutoff=settings['hphob_dist_cut']):
    h1, h2 = interactions.close_contacts(protein.atom_dict[protein.atom_dict['ishydrophobe']], 
                                         ligand.atom_dict[ligand.atom_dict['ishydrophobe']], 
                                         cutoff)

    dist, _type = [], []
    resname_h1, resname_h2 = [], []
    resnum_h1, resnum_h2 = [], []
    atom_h1, atom_h2 = [], []

    for i in range(len(h1)):
        _dist = math.dist(h1['coords'][i], h2['coords'][i])
        dist.append(_dist)
        _resname_h1, _resname_h2 = h1['resname'][i], h2['resname'][i]
        _resnum_h1, _resnum_h2 = h1['resnum'][i], h2['resnum'][i]
        _atom_h1, _atom_h2 = h1['atomtype'][i], h2['atomtype'][i]
        
        if _int == True:
            if _resname_h1 == '':
                _resname_h1 = 'LIG'
            elif _resname_h2 == '':
              _resname_h2 = 'LIG'
        else:
            pass
            
        resname_h1.append(_resname_h1)
        resname_h2.append(_resname_h2)
        
        resnum_h1.append(_resnum_h1)
        resnum_h2.append(_resnum_h2)
        
        atom_h1.append(_atom_h1)
        atom_h2.append(_atom_h2)

        _type.append("hphob")    

    """
    Acceptor - protein ; Donor - ligand
    """    
    
    df = pd.DataFrame({"type"             : _type                          ,
                       "acceptor_resname" : resname_h1                     ,
                       "acceptor_resnum"  : resnum_h1                      ,
                       "acceptor_atom"    : atom_h1                        ,
                       "donor_resname"    : resname_h2                     ,
                       "donor_resnum"     : resnum_h2                      ,
                       "donor_atom"       : atom_h2                        ,
                       "distance"         : dist                           ,
                       "angle"            : [0 for i in range(len(_type))] ,
                       "ring1_resname"    : [0 for i in range(len(_type))] ,
                       "ring1_resnum"     : [0 for i in range(len(_type))] ,
                       "ring2_resname"    : [0 for i in range(len(_type))] ,
                       "ring2_resnum"     : [0 for i in range(len(_type))] ,
                       "cation_atom"      : [0 for i in range(len(_type))] ,
                       })

    return df
    
    
def sbridge_oddt(protein, ligand, _int=True, cutoff=settings['sbridge_dist_cut'], mol1_exact=True, mol2_exact=True):
    p1m2 = interactions.salt_bridge_plus_minus(protein, ligand, cutoff)
    p2m1 = interactions.salt_bridge_plus_minus(ligand, protein, cutoff)
    
    dist_p1m2, dist_p2m1 = [], []
    a_resname_p1m2, c_resname_p1m2 = [], []
    a_resname_p2m1, c_resname_p2m1 = [], []
    a_resnum_p1m2, c_resnum_p1m2 = [], []
    a_resnum_p2m1, c_resnum_p2m1 = [], []
    a_atom_p1m2, c_atom_p1m2 = [], []
    a_atom_p2m1, c_atom_p2m1 = [], []
    _type = []
    
    """
    p1m2 : protein provides cations, ligand provides anions
    """
    for i in range(len(p1m2[0])):
        dist = math.dist(p1m2[0]['coords'][i], p1m2[1]['coords'][i])
        dist_p1m2.append(dist)
        c_resname, a_resname = p1m2[0]['resname'][i], p1m2[1]['resname'][i]
        c_resnum, a_resnum = p1m2[0]['resnum'][i], p1m2[1]['resnum'][i]
        c_atom, a_atom = p1m2[0]['atomtype'][i], p1m2[1]['atomtype'][i]
        
        if _int == True:
            if a_resname == '':
                a_resname = 'LIG'
            elif c_resname == '':
              c_resname = 'LIG'
        else:
            pass
            
        a_resname_p1m2.append(a_resname)
        c_resname_p1m2.append(c_resname)
        
        a_resnum_p1m2.append(a_resnum)
        c_resnum_p1m2.append(c_resnum)
        
        a_atom_p1m2.append(a_atom)
        c_atom_p1m2.append(c_atom)

        _type.append("sbridge")    
    
    """
    p2m1 : protein provides anions, ligand provides cations
    """
    for i in range(len(p2m1[0])):
        dist = math.dist(p2m1[0]['coords'][i], p2m1[1]['coords'][i])
        dist_p2m1.append(dist)
        c_resname, a_resname = p2m1[0]['resname'][i], p2m1[1]['resname'][i]
        c_resnum, a_resnum = p2m1[0]['resnum'][i], p2m1[1]['resnum'][i]
        c_atom, a_atom = p2m1[0]['atomtype'][i], p2m1[1]['atomtype'][i]
        
        if _int == True:
            if a_resname == '':
                a_resname = 'LIG'
            elif c_resname == '':
              c_resname = 'LIG'
        else:
            pass
            
        a_resname_p2m1.append(a_resname)
        c_resname_p2m1.append(c_resname)
        
        a_resnum_p2m1.append(a_resnum)
        c_resnum_p2m1.append(c_resnum)
        
        a_atom_p2m1.append(a_atom)
        c_atom_p2m1.append(c_atom)

        _type.append("sbridge")
      
    """
    Acceptor - protein ; Donor - ligand
    """    
    
    df = pd.DataFrame({"type"             : _type                          ,
                       "acceptor_resname" : c_resname_p1m2 + a_resname_p2m1,
                       "acceptor_resnum"  : c_resnum_p1m2  + a_resnum_p2m1 ,
                       "acceptor_atom"    : c_atom_p1m2    + a_atom_p2m1   ,
                       "donor_resname"    : a_resname_p1m2 + c_resname_p2m1,
                       "donor_resnum"     : a_resnum_p1m2  + c_resnum_p2m1 ,
                       "donor_atom"       : a_atom_p1m2    + c_atom_p2m1   ,
                       "distance"         : dist_p1m2      + dist_p2m1     ,
                       "angle"            : [0 for i in range(len(_type))] ,
                       "ring1_resname"    : [0 for i in range(len(_type))] ,
                       "ring1_resnum"     : [0 for i in range(len(_type))] ,
                       "ring2_resname"    : [0 for i in range(len(_type))] ,
                       "ring2_resnum"     : [0 for i in range(len(_type))] ,
                       "cation_atom"      : [0 for i in range(len(_type))] ,
                       })

    return df
    

def pistack_oddt(protein, ligand, _int=True, cutoff=settings['pistack_dist_cut'], tolerance=30):
    r1, r2 = interactions.close_contacts(protein.ring_dict, ligand.ring_dict, cutoff, x_column='centroid', y_column='centroid')

    dist, _angle, _type = [], [], []
    r1_resname, r1_resnum = [], []
    r2_resname, r2_resnum = [], []
    
    angle1 = angle_2v(r1['vector'], r2['vector'])
    angle2 = angle(r1['vector'] + r1['centroid'],
                   r1['centroid'],
                   r2['centroid'])
    angle3 = angle(r2['vector'] + r2['centroid'],
                   r2['centroid'],
                   r1['centroid'])
        
    for i in range(len(r1)):
        _dist = math.dist(r1[i]['centroid'], r2[i]['centroid'])
        dist.append(_dist)
        _angle.append(angle1[i])
        _type.append("pi_stacking")

        _r1_resnum, _r2_resnum = r1[i]['resnum'], r2[i]['resnum']
        _r1_resname, _r2_resname = r1[i]['resname'], r2[i]['resname']

        if _int == True:
            if _r1_resname == '':
                _r1_resname = 'LIG'
            elif _r2_resname == '':
                _r2_resname = 'LIG'
        else:
            pass
            
        r1_resname.append(_r1_resname)
        r2_resname.append(_r2_resname)
        
        r1_resnum.append(_r1_resnum)
        r2_resnum.append(_r2_resnum)
        
    df = pd.DataFrame({"type"             : _type,
                       "acceptor_resname" : [0 for i in range(len(_type))],
                       "acceptor_resnum"  : [0 for i in range(len(_type))],
                       "acceptor_atom"    : [0 for i in range(len(_type))],
                       "donor_resname"    : [0 for i in range(len(_type))],
                       "donor_resnum"     : [0 for i in range(len(_type))],
                       "donor_atom"       : [0 for i in range(len(_type))],
                       "distance"         : dist,
                       "angle"            : _angle,
                       "ring1_resname"    : r1_resname,
                       "ring1_resnum"     : r1_resnum,
                       "ring2_resname"    : r2_resname,
                       "ring2_resnum"     : r2_resnum,
                       "cation_atom"      : [0 for i in range(len(_type))],
                       })

    return df


def pication_oddt(protein, ligand, _int=True, cutoff=settings['pication_dist_cut'], tolerance=30, cation_exact=False):
    cation_map = ligand.atom_dict['isplus']
    r1, plus2 = interactions.close_contacts(protein.ring_dict, ligand.atom_dict[cation_map], cutoff, x_column='centroid')

    dist, _angle, _type = [], [], []
    r1_resname, r1_resnum = [], []
    plus2_resname, plus2_resnum, plus2_atom = [], [], []
    
    angle1 = angle_2v(r1['vector'], plus2['coords'] - r1['centroid'])

    for i in range(len(r1)):
        _dist = math.dist(r1[i]['centroid'], plus2[i]['coords'])
        dist.append(_dist)
        _angle.append(angle1[i])
        _type.append("pi_cation")

        _r1_resnum, _plus2_resnum = r1[i]['resnum'], plus2[i]['resnum']
        _r1_resname, _plus2_resname = r1[i]['resname'], plus2[i]['resname']
        _plus2_atom = plus2[i]['atomtype']

        if _int == True:
            if _r1_resname == '':
                _r1_resname = 'LIG'
            elif _plus2_resname == '':
                _plus2_resname = 'LIG'
        else:
            pass
            
        r1_resname.append(_r1_resname)
        plus2_resname.append(_plus2_resname)
        
        r1_resnum.append(_r1_resnum)
        plus2_resnum.append(_plus2_resnum)
        
        plus2_atom.append(_plus2_atom)
      

    cation_map = protein.atom_dict['isplus']
    r1, plus2 = interactions.close_contacts(ligand.ring_dict, protein.atom_dict[cation_map], cutoff, x_column='centroid')
    angle1 = angle_2v(r1['vector'], plus2['coords'] - r1['centroid'])

    for i in range(len(r1)):
        _dist = math.dist(r1[i]['centroid'], plus2[i]['coords'])
        dist.append(_dist)
        _angle.append(angle1[i])
        _type.append("pi_cation")

        _r1_resnum, _plus2_resnum = r1[i]['resnum'], plus2[i]['resnum']
        _r1_resname, _plus2_resname = r1[i]['resname'], plus2[i]['resname']

        if _int == True:
            if _r1_resname == '':
                _r1_resname = 'LIG'
            elif _plus2_resname == '':
                _plus2_resname = 'LIG'
        else:
            pass
            
        r1_resname.append(_r1_resname)
        r2_resname.append(_r2_resname)
        
        r1_resnum.append(_r1_resnum)
        r2_resnum.append(_r2_resnum)
        
        plus2_atom.append(_plus2_atom)     
    
    """
    ring2 - cation
    """
    df = pd.DataFrame({"type"             : _type,
                       "acceptor_resname" : [0 for i in range(len(_type))],
                       "acceptor_resnum"  : [0 for i in range(len(_type))],
                       "acceptor_atom"    : [0 for i in range(len(_type))],
                       "donor_resname"    : [0 for i in range(len(_type))],
                       "donor_resnum"     : [0 for i in range(len(_type))],
                       "donor_atom"       : [0 for i in range(len(_type))],
                       "distance"         : dist,
                       "angle"            : _angle,
                       "ring1_resname"    : r1_resname,
                       "ring1_resnum"     : r1_resnum,
                       "ring2_resname"    : plus2_resname,
                       "ring2_resnum"     : plus2_resnum,
                       "cation_atom"      : plus2_atom,
                       })

    return df
    

#######################################################################################
################################### Initialisation ####################################
#######################################################################################

def fingerprint(protein, ligand, _int, settings):
    """
    Parameters
    ----------
    protein, ligand : list, with len(protein) == len(ligand) == 2
        First entry is the Molecule class object based on RDKit.
        Second entry is the oddt.toolkit.Molecule object used exclusively to 
            compute pi_stacking only using ODDT library.
    
    Returns
    -------
    fp : pandas DataFrame, shape = [n rows x 18 columns]
        A table which contains all the computed data, e.g., distance, angles, etc.    
    """
    
    fp1 = hbonds_oddt(protein[0], ligand[0], cutoff=settings['hbond_dist_cut'], _int=True)
    fp2 = hphob_oddt(protein[0], ligand[0], cutoff=settings['hphob_dist_cut'], _int=True)
    fp3 = sbridge_oddt(protein[0], ligand[0], cutoff=settings['sbridge_dist_cut'], _int=True)
    fp4 = pistack_oddt(protein[0], ligand[0], cutoff=settings['pistack_dist_cut'], _int=True)
    fp5 = pication_oddt(protein[0], ligand[0], cutoff=settings['pication_dist_cut'], _int=True)
    fp6 = xbonds_oddt(protein[0], ligand[0], cutoff=settings['xbond_dist_cut'], _int=True)
    
    
    fp = pd.concat([fp1, fp2, fp3, fp4, fp5, fp6], sort=False, ignore_index=True)
    fp.insert(0, "index", int(protein[-1]), True)
    fp.insert(0, "ligand", ligand[-2], True)
    fp.insert(0, "protein", protein[-2], True)
    
    fp['acceptor_resnum'] = fp['acceptor_resnum'].astype(int)
    fp['donor_resnum']    = fp['donor_resnum'].astype(int)
    fp['ring2_resnum']    = fp['ring2_resnum'].astype(int)
    fp['ring1_resnum']    = fp['ring1_resnum'].astype(int)

    
    return fp


#######################################################################################
###################################### User Input #####################################
#######################################################################################

if __name__ == '__main__':
    import sys
    import getopt

    def usage():
        "Print helpful, accurate usage statement to stdout."
        print("Usage: ifp.py -commnand")
        print("")
        print("    Description of command...")
        print("         -p --p  protein files ")
        print("         -l --l  ligand files ")
        print("         -i --i  'pl' for protein-ligand  interaction")
        print("                 'pp' for protein-protein interaction")
        print("         -h --h  get_help ")

    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'p:l:i:h')

    except getopt.GetoptError:
        print("Arguments Error!")
        usage()
        sys.exit(2)

    # initialize required parameters
    # -p: protein
    # -l: ligand
    # -i: interaction
    protein_name = None
    ligand_name =  None
    interaction = None

    #'p:l:i:h'
    for o, a in opt_list:
        if o in ('-p', '--p'):
            protein_name = str(a)
        if o in ('-l', '--l'):
            ligand_name = str(a)
        if o in ('-i', '--i'):
            interaction = str(a)
        if o in ('-h', '--h'):
            usage()
            sys.exit()


    if not protein_name:
        print('ifp.py: Protein file must be specified.')
        usage()
        sys.exit()
        
    if not ligand_name:
        print('ifp.py: Ligand file must be specified.')
        usage()
        sys.exit()
        
    if not interaction:
        print('ifp.py: Interaction type must be specified.')
        usage()
        sys.exit()
     
     
    if interaction not in ['pl', 'pp']:
        print("Interaction type error!")  
        print("Must be either 'pl' or 'pp'! ")    
        usage()
        sys.exit()
        
         
    #######################################################################################
    ################################### Main function #####################################
    #######################################################################################
    
    print(" ")
    print("Initialising ...")
    fps = pd.DataFrame()

    print("*******************************") 

    path = r"./"
    protein_file = path + "protein/" + protein_name 
    ligand_file = path + "ligand/" + ligand_name 
    pro_name = protein_file.split("/")[-1].split(".")[0]
    lig_name = ligand_file.split("/")[-1].split(".")[0]
    
    print("Reading protein file {}".format(pro_name))
    print("Reading ligand file {}".format(lig_name))
    
    protein_oddt = next(oddt.toolkit.readfile('pdb', protein_file))
    protein_oddt.protein=True
    protein_oddt.addh()
                
    mols = list(oddt.toolkit.readfile('pdb', ligand_file))
    ligand_oddt = mols[0]
    ligand_oddt.addh()
    
    if interaction == 'pl':
        _int = True
    elif interaction == 'pp':
        _int = False
    else:
        print("Interaction type error! Refer to the usage panel.")
            
    protein_duo = [protein_oddt, pro_name, 0]
    ligand_duo  = [ligand_oddt , lig_name, 0]
            
    tmp  = fingerprint(protein_duo, ligand_duo, _int, settings)
    fps = pd.concat([fps, tmp], sort=False)

    print(" ")
    print("Writing descriptors file ...")
    
    #fps = fps[fps["protein_res"].str.contains("HOH") == False]
    #fps = fps[fps['dist'] <= 6.5] 

    print(fps.head(10))

    fps.to_csv('./descriptors/descriptors.csv', index=False)
    
    print("Process complete.")

