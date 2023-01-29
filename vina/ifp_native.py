import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import oddt
from oddt import interactions
from scipy.spatial.distance import cdist

from rdkit import Chem
from rdkit.Chem.rdmolfiles import MolFromPDBFile
from rdkit.Chem.rdmolfiles import MolFromSmarts
from rdkit.Chem.rdmolops import AddHs

import numpy as np
import pandas as pd

import os
import glob
import math



################################################################################
################################# Fundamentals #################################
################################################################################

def resname(atom):
    info = atom.GetPDBResidueInfo()
    if info is None:
        return ''
    return ':'.join(map(lambda x: str(x).strip(),
                        [info.GetChainId(), str(info.GetResidueNumber()), info.GetResidueName(), info.GetInsertionCode()]))

def atomname(atom):
    pdb = atom.GetPDBResidueInfo()
    if pdb is None:
        return str(atom.GetIdx())
    return pdb.GetName().strip()

def coords(atom):
    return atom.GetOwningMol().GetConformer(0).GetAtomPosition(atom.GetIdx())

def distance(atom1, atom2):
    return coords(atom1).Distance(coords(atom2))
    
def distance_oddt(x, y):
    """Computes distance between each pair of points from x and y.
       Used exclusively for pi stacking from ODDT library only.
    Parameters
    ----------
    x : numpy arrays, shape = [n_x, 3]
        Array of poinds in 3D
    y : numpy arrays, shape = [n_y, 3]
        Array of poinds in 3D
    
    Returns
    -------
    dist_matrix : numpy arrays, shape = [n_x, n_y]
        Distance matrix
    """
    return cdist(x, y)
    
def centroid_coords(atoms):
    _coords = np.array([coords(atom) for atom in atoms])
    _coords = _coords.mean(axis=0)
    return _coords

def angle_atom(atom1, atom2, atom3):
    v1 = coords(atom1) - coords(atom2)
    v3 = coords(atom3) - coords(atom2)
    return v1.AngleTo(v3) * 180.0 / np.pi

def angle_vector(v1, v2):
    v1 /= np.linalg.norm(v1)
    v2 /= np.linalg.norm(v2)
    angle =  np.arccos(np.clip(np.dot(v1, v2), -1.0, 1.0))
    angle *= 180 / np.pi
    if angle > 90:
        angle = 180 - angle
    assert 0 <= angle <= 90, angle
    return angle

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
    
def _check_angles(angles, hybridizations, tolerance):
    """Helper function for checking if interactions are strict"""
    angles = np.nan_to_num(angles)  # NaN's throw warning on comparisons
    ideal_angles = np.take(BASE_ANGLES, hybridizations)[:, np.newaxis]
    lower_bound = ideal_angles - tolerance
    upper_bound = ideal_angles + tolerance
    return ((angles > lower_bound) & (angles < upper_bound)).any(axis=-1)
        
class Molecule:
    """Initialise various information needed for the computation and 
       return as attributes of the input molecule.
    Parameters
    ----------
    mol : RDKit Mol object
        Either protein or ligand read by MolFromPDBFile(), explicit hydrogens added.
    is_protein : bool
        Is mol protein, or not
    settings : dict
        Carries information such as cutoff distance and angles for interactions.
        
    Returns
    -------
    Vary for different interactions.   
    """
    
    def __init__(self, mol, is_protein, settings):
        self.mol = mol
        self.is_protein = is_protein
        self.settings = settings

        self.hbond_donors, self.hbond_acceptors = self.init_hbond()
        self.charged, self.charge_groups = self.init_saltbridge()
        self.contacts = self.init_contacts()

    ################################# Hydrogen bonds ##################################
    
    def init_hbond(self):
        """ 
        Initialise hydrogen bond donor and acceptor atoms, `donors` and `acceptors`,
        respectively, via _is_donor() and _is_acceptor() functions, respectively,
        as defined below.
        
        Note that we only consider nitrogen and oxygen as hydrogen bond donors and
        acceptors.
        
        See criteria below.
        """    
        donors = [atom for atom in self.mol.GetAtoms() if self._is_donor(atom)]
        acceptors = [atom for atom in self.mol.GetAtoms() if self._is_acceptor(atom)]
        return donors, acceptors

    def _is_donor(self, atom):
        if atom.GetAtomicNum() in [7, 8]:
            if _get_bonded_hydrogens(atom):
                return True
        return False

    def _is_acceptor(self, atom):
        if atom.GetAtomicNum() == 8:
            return True
        if atom.GetAtomicNum() == 7 and atom.GetExplicitValence() < 4:
            return True
        return False

    ################################### Salt bridges ##################################        
    
    def init_saltbridge(self):
        """ 
        Initialise charged atoms `charged` and charged residues `charge_groups` of
        both protein and ligand via _symmetric_charged_protein_atoms() and
        symmetric_charged_ligand_atoms(), respectively, as defined below.
        """    
        charged = [atom for atom in self.mol.GetAtoms() if atom.GetFormalCharge() != 0]
        
        if self.is_protein:
            charge_groups = self._symmetric_charged_protein_atoms()
            
        else:
            charge_groups = self._symmetric_charged_ligand_atoms()
            
        return charged, charge_groups
        
    def _symmetric_charged_protein_atoms(self):
        protein_groups = {}
        
        for protein_atom in self.mol.GetAtoms():
            if atomname(protein_atom) in ["OD1", "OD2", "OE1", "OE2", "NH1", "NH2"]:
                if resname(protein_atom) not in protein_groups:
                    protein_groups[resname(protein_atom)] = []
                protein_groups[(resname(protein_atom))] += [protein_atom]
        
        return protein_groups
        
    def _symmetric_charged_ligand_atoms(self):
        ligand_groups = {}
        smartss = [('[CX3](=O)[O-]', 2, [1, 2]),
                   ('[CX3](=[NH2X3+])[NH2X3]', 1, [1, 2])]
        
        idx_to_atom = {atom.GetIdx(): atom for atom in self.mol.GetAtoms()}
        
        for smarts, k, v in smartss:
            mol = MolFromSmarts(smarts)
            matches = self.mol.GetSubstructMatches(mol)
            
            for match in matches:
                ligand_groups[match[k]] = [idx_to_atom[match[_v]] for _v in v]
                
        return ligand_groups
        
    ############################## Hydrophobic contacts ###############################
    
    def init_contacts(self):
        """
        Initialise coordinates `coord`, van der Waals radius `vdw`,
        residue names `res_name` and atom names `atom_name` for 
        computing hydrophobic contacts and atom id `atom_id` for
        building bitstring.
        """
        coord, vdw, atom_name, res_name, atom_id = [], [], [], [], []
        
        for atom in self.mol.GetAtoms():
            atom_id += [atom.GetIdx()]
            
            if atom.GetAtomicNum() not in self.settings['nonpolar']:
                continue
            
            coord += [coords(atom)]
            atom_name += [atomname(atom)]
            res_name += [resname(atom)]
            vdw += [self.settings['nonpolar'][atom.GetAtomicNum()]]
            
        if coord:
            coord = np.vstack(coord)
            vdw = np.array(vdw)
        
        return coord, vdw, res_name, atom_name, atom_id
        
    def get_aromatic_rings(self):
        return [ring for ring in self.mol.GetRingInfo().AtomRings()
                if self.mol.GetAtomWithIdx(ring[0]).GetIsAromatic()]
                
    def get_centroid(self, atom_idx):
        atoms = [self.mol.GetAtomWithIdx(a) for a in atom_idx]
        return centroid_coords(atoms)
        
    def get_normal(self, ring):
        centroid = self.get_centroid(ring)
        coords1 = coords(self.mol.GetAtomWithIdx(ring[0])) - centroid
        coords2 = coords(self.mol.GetAtomWithIdx(ring[1])) - centroid
        
        normal = np.cross(coords1, coords2)
        normal /= np.linalg.norm(normal)
        
        return normal    
        
              
#######################################################################################
########################### Compute atom-level interactions ###########################
#######################################################################################

def _get_bonded_hydrogens(atom):
    hydrogens = []
    for bond in atom.GetBonds():
        if bond.GetBeginAtomIdx() != atom.GetIdx():
            hydrogen = bond.GetBeginAtom()
        else:
            hydrogen = bond.GetEndAtom()
            
        if hydrogen.GetAtomicNum() == 1:
            hydrogens += [hydrogen]
    return hydrogens


def _hbond_hydrogen_angle(acceptor, donor):
    best_angle, best_hydrogen = 0, None
    for hydrogen in _get_bonded_hydrogens(donor):
        _angle = angle_atom(donor, hydrogen, acceptor)
        if _angle > best_angle:
            best_angle = _angle
            best_hydrogen = hydrogen
    return best_hydrogen, best_angle


def _hbond_compute(donor_mol, acceptor_mol, settings, protein_is_donor, bitstring):
    """
    Parameters
    ----------
    donor_mol, acceptor_mol : Molecule class objects
        Hydrogen donor and acceptor molecules based on RDKit Mol.
    settings: dict
        Contain cutoff distances and angles.
    protein_is_donor: bool
        Whether protein is always the hydrogen bonds donor
    
    Returns
    -------
    hbonds : pandas DataFrame
        Information related to the computed hydrogen bonds.
    """
    hbonds = []
    
    for donor in donor_mol.hbond_donors:
        for acceptor in acceptor_mol.hbond_acceptors:
            for hydrogen in _get_bonded_hydrogens(donor):
                dist = distance(acceptor, hydrogen)
                
                if dist > settings['hbond_dist_cut']: 
                    continue
                
                angle = angle_atom(donor, hydrogen, acceptor)
                
                if angle < settings['hbond_angle_cut']: 
                    continue

                if protein_is_donor:
                    label = 'hbond_donor'
                    protein_atom = donor
                    ligand_atom = acceptor

                    bitstring[donor.GetIdx()] = bitstring[donor.GetIdx()][0:2] + '1' + bitstring[donor.GetIdx()][3:]
                    
                else:
                    label = 'hbond_acceptor'
                    protein_atom = acceptor
                    ligand_atom = donor
                    
                    bitstring[acceptor.GetIdx()] = bitstring[acceptor.GetIdx()][0:1] + '1' + bitstring[acceptor.GetIdx()][2:]
                    

                hbonds += [{'label': label,
                            'protein_res': resname(protein_atom),
                            'protein_atom': atomname(protein_atom),
                            'ligand_atom': atomname(ligand_atom),
                            'dist': dist,
                            'angle': angle,
                            'hydrogen': atomname(hydrogen)}]
    
    return hbonds, bitstring


def hbond_compute(protein, ligand, settings, bitstring):
    """
    Parameters
    ----------
    protein, ligand : Molecule class objects
        Protein (pocket) and output ligand based on RDKit Mol.
    settings: dict
        Contain cutoff distances and angles.
    
    Returns
    -------
    acceptor + donor : pandas DataFrame
        Information related to the hydrogen bond acceptors and donors.
    """
    
    
    donor, bitstring = _hbond_compute(protein, ligand, settings, True, bitstring)
    acceptor, bitstring = _hbond_compute(ligand, protein, settings, False, bitstring)
    
    return acceptor + donor, bitstring
    
    
def saltbridge_compute(protein, ligand, settings, bitstring):
    """
    Parameters
    ----------
    protein, ligand : Molecule class objects
        Protein (pocket) and output ligand based on RDKit Mol.
    settings: dict
        Contain cutoff distances.
    
    Returns
    -------
    contacts : pandas DataFrame
        Information related to the computed salt bridges.
    """
    
    saltbridges = []
    
    """
    Original authors note
    ---------------------
    Note that much of the complexity here stems from taking into account
    symetric atoms. Specifically for carboxylate and guanidinium groups,
    we consider not just the atom that is arbitrarily assigned a formal
    charge, but also the atom that is charged in the other resonance
    structure.
    """
 
    for protein_atom in protein.charged:
        for ligand_atom in ligand.charged:
            lig_charge = ligand_atom.GetFormalCharge()
            protein_charge = protein_atom.GetFormalCharge()
            if lig_charge * protein_charge >= 0: continue

            # Expand protein_atom and ligand_atom to all symetric atoms
            # ... think carboxylates and guanidiniums.
            if ('saltbridge_resonance' in settings and ligand_atom.GetIdx() in ligand.charge_groups):
                ligand_atoms = ligand.charge_groups[ligand_atom.GetIdx()]
            else:
                ligand_atoms = [ligand_atom]
            
            if ('saltbridge_resonance' in settings and resname(protein_atom) in protein.charge_groups):
                protein_atoms = protein.charge_groups[resname(protein_atom)]
            else:
                protein_atoms = [protein_atom]

            # Get minimum distance between any pair of protein and ligand
            # atoms in the groups.
            dist = float('inf')
            for _ligand_atom in ligand_atoms:
                for _protein_atom in protein_atoms:
                    _dist = distance(_protein_atom, _ligand_atom)
                    if _dist < dist:
                        dist = _dist
                        closest_protein_atom = _protein_atom
                        closest_ligand_atom = _ligand_atom
            
            bitstring[_protein_atom.GetIdx()] = bitstring[_protein_atom.GetIdx()][0:3] + '1' + bitstring[_protein_atom.GetIdx()][4:]
            bitstring[_ligand_atom.GetIdx()] = bitstring[_ligand_atom.GetIdx()][0:3] + '1' + bitstring[_ligand_atom.GetIdx()][4:]
            
            if dist < settings['sbridge_dist_cut']:
                saltbridges += [{'label': 'saltbridge',
                                 'protein_res': resname(closest_protein_atom),
                                 'protein_atom': atomname(closest_protein_atom),
                                 'ligand_atom': atomname(closest_ligand_atom),
                                 'dist': dist}]
    
    return saltbridges, bitstring      


def contact_compute(protein, ligand, settings, bitstring):
    """
    Parameters
    ----------
    protein, ligand : Molecule class object
        Protein (pocket) and output ligand based on RDKit Mol.
    settings: dict
        Contain cutoff distances.
    
    Returns
    -------
    contacts : pandas DataFrame
        Information related to the computed hydrophobic contacts.
    """
    protein = protein.contacts
    ligand = ligand.contacts
    
    dists = protein[0].reshape(1, -1, 3) - ligand[0].reshape(-1, 1, 3)
    dists = np.linalg.norm(dists, axis=2)
    vdw = protein[1].reshape(1, -1) + ligand[1].reshape(-1, 1)
    contact_idx = np.argwhere(dists < vdw * settings['contact_scale_cut'])
    
    contacts = []

    for i, j in contact_idx:
        contacts += [{'label': 'contact',
                      'protein_res': protein[2][j],
                      'protein_atom': protein[3][j],
                      'ligand_atom': ligand[3][i],
                      'dist': dists[i, j],
                      'vdw': vdw[i, j]}]
                      
        bitstring[int(protein[4][j])] = bitstring[int(protein[4][j])][0:4] + '1' + bitstring[int(protein[4][j])][5:]
                      
    return contacts, bitstring


def close_contacts(x, y, cutoff, x_column='coords', y_column='coords', cutoff_low=0.):
    """Returns pairs of atoms which are within close contac distance cutoff.
    The cutoff is semi-inclusive, i.e (cutoff_low, cutoff].
    Parameters
    ----------
    x, y : atom_dict-type numpy array
        Atom dictionaries generated by oddt.toolkit.Molecule objects.
    cutoff : float
        Cutoff distance for close contacts
    x_column, y_column : string, (default='coords')
        Column containing coordinates of atoms (or pseudo-atoms,
        i.e. ring centroids)
    cutoff_low : float (default=0.)
        Lower bound of contacts to find (exclusive). Zero by default.
        .. versionadded:: 0.6
    
    Returns
    -------
    x_, y_ : atom_dict-type numpy array
        Aligned pairs of atoms in close contact for further processing.
    """
    if len(x[x_column]) > 0 and len(x[x_column]) > 0:
        d = distance_oddt(x[x_column], y[y_column])
        index = np.argwhere((d > cutoff_low) & (d <= cutoff))
        return x[index[:, 0]], y[index[:, 1]]
    
    else:
        return x[[]], y[[]]

def _pistack_compute(mol1, mol2, cutoff, tolerance):
    """Returns pairs of rings, which meet pi stacking criteria
    Parameters
    ----------
    mol1, mol2 : oddt.toolkit.Molecule object
        Molecules to compute ring pairs
    cutoff : float, (default=5)
        Distance cutoff for Pi-stacking pairs
    tolerance : int, (default=30)
        Range (+/- tolerance) from perfect direction (parallel or
        perpendicular) in which pi-stackings are considered as strict.
    
    Returns
    -------
    r1, r2 : ring_dict-type numpy array
        Aligned arrays of rings forming pi-stacking
    strict_parallel : numpy array, dtype=bool
        Boolean array align with ring pairs, informing whether rings
        form 'strict' parallel pi-stacking. If false, only distance cutoff is met,
        therefore the stacking is 'crude'.
    strict_perpendicular : numpy array, dtype=bool
        Boolean array align with ring pairs, informing whether rings
        form 'strict' perpendicular pi-stacking (T-shaped, T-face, etc.).
        If false, only distance cutoff is met, therefore the stacking is 'crude'.
    """
    r1, r2 = close_contacts(mol1.ring_dict, mol2.ring_dict, cutoff,
                            x_column='centroid', y_column='centroid')
    
    if len(r1) > 0 and len(r2) > 0:
        angle1 = angle_2v(r1['vector'], r2['vector'])
        angle2 = angle(r1['vector'] + r1['centroid'],
                       r1['centroid'],
                       r2['centroid'])
        angle3 = angle(r2['vector'] + r2['centroid'],
                       r2['centroid'],
                       r1['centroid'])
        
        strict_parallel = (((angle1 > 180 - tolerance) | (angle1 < tolerance)) &
                           ((angle2 > 180 - tolerance) | (angle2 < tolerance) |
                            (angle3 > 180 - tolerance) | (angle3 < tolerance)))
        strict_perpendicular = ((angle1 > 90 - tolerance) & (angle1 < 90 + tolerance) &
                                (((angle2 > 180 - tolerance) | (angle2 < tolerance)) &
                                 ((angle3 > 90 - tolerance) | (angle3 < 90 + tolerance)) |
                                 ((angle2 > 90 - tolerance) | (angle2 < 90 + tolerance)) &
                                 ((angle3 > 180 - tolerance) | (angle3 < tolerance))
                                 )
                                )
        return r1, r2, strict_parallel, strict_perpendicular, angle1, angle2, angle3
    
    else:
        return r1, r2, np.array([], dtype=bool), np.array([], dtype=bool), np.array([], dtype=bool), np.array([], dtype=bool), np.array([], dtype=bool)
        
            
# Finding interactions using `oddt` library        
def pistack_compute(protein, ligand, settings, bitstring):
    """Returns DataFrame that contains all the computed pi stacking interactions
    Parameters
    ----------
    protein, ligand : oddt.toolkit.Molecule object
        Molecules to compute pi stacking interactions
    settings : dict
        Contains two important information - `cutoff` and `tolerance`:
        
        -> cutoff : float, (default=5)
               Distance cutoff for Pi-stacking pairs.
        -> tolerance : int, (default=30)
               Range (+/- tolerance) from perfect direction (parallel or
               perpendicular) in which pi-stackings are considered as strict.
    
    Returns
    -------
    pistack : pandas DataFrame
        Information related to the computed pi stacking interactions.
    """    
    
    pistack, atom_id = [], []
    
    ring1, ring2, strict_para, strict_perp, angle1, angle2, angle3 = _pistack_compute(protein, ligand, 
                                                                                      cutoff=settings['pistack_dist_cut'], 
                                                                                      tolerance=settings['pistack_angle_cut'])
                                                                           
    """
    ring1, ring2 : ring_dict-type numpy array
        Aligned arrays of rings forming pi-stacking
    strict_parallel : numpy array, dtype=bool
        Boolean array align with ring pairs, informing whether rings
        form 'strict' parallel pi-stacking. If false, only distance cutoff is met,
        therefore the stacking is 'crude'.
    strict_perpendicular : numpy array, dtype=bool
        Boolean array align with ring pairs, informing whether rings
        form 'strict' perpendicular pi-stacking (T-shaped, T-face, etc.).
        If false, only distance cutoff is met, therefore the stacking is 'crude'.
    angle1, angle2, angle3 : float
        angle1 - the angle between two vectors which are perpendicular to the ring plane
        angle2 - the angle between (ring1_member_atom, ring1_centroid, ring2_centroid)
        angle3 - the angle between (ring2_member_atom, ring2_centroid, ring1_centroid)
    """
    
    ring1_resnum, ring1_resname = ring1['resnum'].tolist(), ring1['resname'].tolist()
        
    ring1_centroid, ring2_centroid = ring1['centroid'].tolist(), ring2['centroid'].tolist()
     
    assert(len(ring1_centroid) == len(ring2_centroid))    
    r1, r2 = close_contacts(protein.ring_dict, ligand.ring_dict, 
                            cutoff=settings['pistack_dist_cut'], 
                            x_column='centroid', y_column='centroid')
                            
    for j in range(len(ring1['resid'])):
        num_atoms = len(protein.residues[int(ring1['resid'][j])].atoms)
        
        for k in range(num_atoms):
            atom_id += [protein.residues[int(ring1['resid'][j])].atoms[k].idx0]
       
    for _atom_id in atom_id:
        bitstring[_atom_id] = bitstring[_atom_id][0:5] + '1' + bitstring[_atom_id][6:]
        
    
    for i in range(len(ring1_centroid)):
        r1_centroid = ring1['centroid'][i]
        r2_centroid = ring2['centroid'][i]
        
        dist = math.dist(r1_centroid, r2_centroid)
        
        _angle1, _angle2, _angle3 = angle1[i], angle2[i], angle3[i]
            
        r1_res = ":" + str(ring1_resnum[i]) + ":" + str(ring1_resname[i]) + ":"
        
        isPara, isPerp = strict_para[i], strict_perp[i]
        
        #bitstring[ring1['id'][i]] = bitstring[ring1['id'][i]][0:5] + '1' + bitstring[ring1['id'][i]][6:]
        bitstring = bitstring
    
        pistack += [{'label': 'pi_stacking',
                     'protein_res': r1_res,
                     'protein_atom': '',
                     'ligand_atom': '',
                     'dist': dist,
                     'angle': round(float(_angle1),2),
                     'hydrogen':'',
                     'vdw': '',
                     'protein_ring_centroid_x': r1_centroid[0],
                     'protein_ring_centroid_y': r1_centroid[1],
                     'protein_ring_centroid_z': r1_centroid[2],
                     'ligand_ring_centroid_x': r2_centroid[0],
                     'ligand_ring_centroid_y': r2_centroid[1],
                     'ligand_ring_centroid_z': r2_centroid[2],
                     'angle_ctr1': round(float(_angle2),2),
                     'angle_ctr2': round(float(_angle3),2),
                     'is_Parallel': isPara,
                     'is_Perpendicular': isPerp}]
        
        
    return pistack, bitstring
    
def _pication_compute(mol1, mol2, cutoff=5, tolerance=30, cation_exact=False):
    """Returns pairs of ring-cation atoms, which meet pi-cation criteria

    Parameters
    ----------
    mol1, mol2 : oddt.toolkit.Molecule object
        Molecules to compute ring-cation pairs

    cutoff : float, (default=5)
        Distance cutoff for Pi-cation pairs

    tolerance : int, (default=30)
        Range (+/- tolerance) from perfect direction (perpendicular)
        in which pi-cation are considered as strict.

    cation_exact : bool
        Requires interacting atoms to have non-zero formal charge.

    Returns
    -------
    r1 : ring_dict-type numpy array
        Aligned rings forming pi-stacking

    plus2 : atom_dict-type numpy array
        Aligned cations forming pi-cation

    strict_parallel : numpy array, dtype=bool
        Boolean array align with ring-cation pairs, informing whether
        they form 'strict' pi-cation. If false, only distance cutoff is met,
        therefore the interaction is 'crude'.

    """
    cation_map = mol2.atom_dict['isplus']
    
    if cation_exact:
        cation_map = cation_map & (mol2.atom_dict['formalcharge'] > 0)
    
    r1, plus2 = close_contacts(mol1.ring_dict,
                               mol2.atom_dict[cation_map],
                               cutoff,
                               x_column='centroid')
    
    if len(r1) > 0 and len(plus2) > 0:
        angle1 = angle_2v(r1['vector'], plus2['coords'] - r1['centroid'])
        ideal_angle = 30  # angle to normal vector
        strict = (
            ((angle1 > ideal_angle - tolerance) & (angle1 < ideal_angle + tolerance)) |
            ((angle1 > 180 - ideal_angle - tolerance) & (angle1 < 180 - ideal_angle + tolerance))
        )
        return r1, plus2, strict, angle1
    
    else:
        return r1, plus2, np.array([], dtype=bool), np.array([], dtype=bool)

       
def pication_compute(protein, ligand, settings, bitstring):
    """Returns DataFrame that contains all the computed pi stacking interactions
    Parameters
    ----------
    protein, ligand : oddt.toolkit.Molecule object
        Molecules to compute pi stacking interactions
    settings : dict
        Contains two important information - `cutoff` and `tolerance`:
        
        -> cutoff : float, (default=5)
               Distance cutoff for Pi-stacking pairs.
        -> tolerance : int, (default=30)
               Range (+/- tolerance) from perfect direction (parallel or
               perpendicular) in which pi-stackings are considered as strict.
    
    Returns
    -------
    pication : pandas DataFrame
        Information related to the computed pi-cation interactions.
    """    
    
    pication, atom_id = [], []
    
    ring, cation, strict, angle = _pication_compute(protein, 
                                                    ligand, 
                                                    cutoff=settings['pication_dist_cut'], 
                                                    tolerance=settings['pication_angle_cut'])
                                                                           
    """
    ring, cation : ring_dict-type numpy array
        Aligned arrays of rings and cations forming pi-cation
    strict : numpy array, dtype=bool
        Boolean array align with ring-cation pairs, informing whether
        they form 'strict' pi-cation. If false, only distance cutoff is met,
        therefore the interaction is 'crude'.
    angle : float
        angle1 - the angle between two vectors which are perpendicular to the ring plane
        angle2 - the angle between (ring1_member_atom, ring1_centroid, ring2_centroid)
        angle3 - the angle between (ring2_member_atom, ring2_centroid, ring1_centroid)
    """
    
    for j in range(len(ring['resid'])):
        num_atoms = len(protein.residues[int(ring['resid'][j])].atoms)
        
        for k in range(num_atoms):
            atom_id += [protein.residues[int(ring['resid'][j])].atoms[k].idx0]
    
    
    for _atom_id in atom_id:
        bitstring[_atom_id] = bitstring[_atom_id][0:6] + '1' + bitstring[_atom_id][7:]
        
    
    for i in range(len(ring["centroid"])):
        ring_centroid = ring["centroid"][i]
        cation_coords = cation["coords"][i]
        
        
        dist = math.dist(ring_centroid, cation_coords)
        _angle = angle[i]
        
        if ring["resname"][i] == '':
            r_res = ":" + str(ring["resnum"][i]) + ":" + "LIG" + ":"
        else:
            r_res = ":" + str(ring["resnum"][i]) + ":" + str(ring["resname"][i]) + ":"
        
        isStrict = strict[i]
    
        pication += [{'label': 'pi_cation',
                      'protein_res': r_res,
                      'protein_atom': '',
                      'ligand_atom': '',
                      'dist': dist,
                      'angle': round(float(_angle),2),
                      'hydrogen':'',
                      'vdw': cation[i]["radius"],
                      'protein_ring_centroid_x': ring_centroid[0],
                      'protein_ring_centroid_y': ring_centroid[1],
                      'protein_ring_centroid_z': ring_centroid[2],
                      'ligand_ring_centroid_x': cation_coords[0],
                      'ligand_ring_centroid_y': cation_coords[1],
                      'ligand_ring_centroid_z': cation_coords[2],
                      'angle_ctr1': '',
                      'angle_ctr2': '',
                      'is_Parallel': '',
                      'is_Perpendicular': '',
                      'charge': cation[i]["charge"],
                      'cation_atom': cation[i]["atomicnum"],
                      }]
        
    return pication, bitstring


def halogenbond_acceptor_halogen(mol1,
                                 mol2,
                                 tolerance=30,
                                 cutoff=4):
    """Returns pairs of acceptor-halogen atoms, which meet halogen bond criteria
    Parameters
    ----------
    mol1, mol2 : oddt.toolkit.Molecule object
        Molecules to compute halogen bond acceptor and halogen pairs
    cutoff : float, (default=4)
        Distance cutoff for A-H pairs
    tolerance : int, (default=30)
        Range (+/- tolerance) from perfect direction defined by atoms hybridization
        in which halogen bonds are considered as strict.
    Returns
    -------
    a, h : atom_dict-type numpy array
        Aligned arrays of atoms forming halogen bond, firstly acceptors,
        secondly halogens
    strict : numpy array, dtype=bool
        Boolean array align with atom pairs, informing whether atoms
        form 'strict' halogen bond (pass all angular cutoffs). If false,
        only distance cutoff is met, therefore the bond is 'crude'.
    """
    a, h = close_contacts(mol1.atom_dict[mol1.atom_dict['isacceptor']],
                          mol2.atom_dict[mol2.atom_dict['ishalogen']],
                          cutoff)
                          
    # skip empty values
    if len(a) > 0 and len(h) > 0:
        angle1 = angle(h['coords'][:, np.newaxis, :],
                       a['coords'][:, np.newaxis, :],
                       a['neighbors'])
        angle2 = angle(a['coords'][:, np.newaxis, :],
                       h['coords'][:, np.newaxis, :],
                       h['neighbors'])
        strict = (_check_angles(angle1, a['hybridization'], tolerance) &
                  _check_angles(angle2, np.ones_like(h['hybridization']), tolerance))
        
        angle_haan, angle_ahhn = [], []
        for i in range(len(h)):
            angle_ahhn.append(angle2[i][0])
            angle_haan.append(angle1[i][0])
        
        return a, h, strict, np.array(angle_haan), np.array(angle_ahhn)
    else:
        return a, h, np.array([], dtype=bool), np.array([], dtype=bool), np.array([], dtype=bool)


def _Xbond_compute(mol1, mol2, cutoff=4, tolerance=30):
    """Calculates halogen bonds between molecules
    Parameters
    ----------
    mol1, mol2 : oddt.toolkit.Molecule object
        Molecules to compute halogen bond acceptor and halogen pairs
    cutoff : float, (default=4)
        Distance cutoff for A-H pairs
    tolerance : int, (default=30)
        Range (+/- tolerance) from perfect direction defined by atoms hybridization
        in which halogen bonds are considered as strict.
    Returns
    -------
    mol1_atoms, mol2_atoms : atom_dict-type numpy array
        Aligned arrays of atoms forming halogen bond
    strict : numpy array, dtype=bool
        Boolean array align with atom pairs, informing whether atoms
        form 'strict' halogen bond (pass all angular cutoffs). If false,
        only distance cutoff is met, therefore the bond is 'crude'.
    """
    a1, h1, s1, g11, g12 = halogenbond_acceptor_halogen(mol1, mol2, cutoff=cutoff, tolerance=tolerance)
    a2, h2, s2, g21, g22 = halogenbond_acceptor_halogen(mol2, mol1, cutoff=cutoff, tolerance=tolerance)
    
    return np.concatenate((a1, h2)), np.concatenate((h1, a2)), np.concatenate((s1, s2)), np.concatenate((g11.reshape(-1,1), g21.reshape(-1,1))), np.concatenate((g12.reshape(-1,1), g22.reshape(-1,1)))
    
    
def Xbond_compute(protein, ligand, settings, bitstring):
    
    Xbond, atom_id = [], []
    
    mol1, mol2, strict, _angle_haan, _angle_ahhn = _Xbond_compute(protein, 
                                                                  ligand, 
                                                                  cutoff=settings['Xbond_dist_cut'], 
                                                                  tolerance=settings['Xbond_angle_cut'])
                                                  
    for j in range(len(mol1['resid'])):
        num_atoms = len(protein.residues[int(mol1['resid'][j])].atoms)
        
        for k in range(num_atoms):
            atom_id += [protein.residues[int(mol1['resid'][j])].atoms[k].idx0]
    
    
    for _atom_id in atom_id:
        bitstring[_atom_id] = bitstring[_atom_id][0:7] + '1'
    
    
    for i in range(len(mol1)):
        mol1_coords = mol1['coords'][i]
        mol2_coords = mol2['coords'][i]
        
        dist = math.dist(mol1_coords, mol2_coords)
        angle_haan, angle_ahhn = _angle_haan[i], _angle_ahhn[i]
        
        mol1_res = ":" + str(mol1["resnum"][i]) + ":" + str(mol1["resname"][i]) + ":"
        #mol2_res = ":" + str(mol2["resnum"][i]) + ":" + "LIG" + ":"        
        
        isStrict = strict[i]
        
        #bitstring[mol1['id'][i]] = bitstring[mol1['id'][i]][0:7] + '1' 

        Xbond += [{'label': 'Xbond',
                   'protein_res': mol1_res,
                   'protein_atom': mol1["atomtype"][i],
                   'ligand_atom': mol2["atomtype"][i],
                   'dist': dist,
                   'angle':'',
                   'hydrogen':'',
                   'vdw': mol2["radius"][i],
                   'protein_ring_centroid_x': '',
                   'protein_ring_centroid_y': '',
                   'protein_ring_centroid_z': '',
                   'ligand_ring_centroid_x': '',
                   'ligand_ring_centroid_y': '',
                   'ligand_ring_centroid_z': '',
                   'angle_ctr1': angle_haan[0],
                   'angle_ctr2': angle_ahhn[0],
                   'is_Parallel': '',
                   'is_Perpendicular': '',
                   'charge': mol2["charge"][i],
                   'cation_atom': '',
                   'halogen_atom': mol2["atomicnum"][i],
                   }]

    return Xbond, bitstring

   
#######################################################################################
############################ Compute residue-level scores #############################
#######################################################################################

"""
def _piecewise_hbond(data, bound=1.2, angle_cut=90.0, mu=2.35, sigma=0.69):

    hbond_acceptor = data[data["label"].str.contains("hbond_acceptor") == True]
    hbond_donor = data[data["label"].str.contains("hbond_donor") == True]
    
    
    def _fdist(x, mu, sigma):
        return (1/(sigma*np.sqrt(2*np.pi)))*(np.exp(-0.5*((x-mu)/sigma)**2))
    
    def fdist(x, bound, mu, sigma):
        return (_fdist(x, mu, sigma)-_fdist(bound, mu, sigma))/(_fdist(mu, mu, sigma))
    
    def fangle(x, angle_cut):
        if (90.0 <= x) and (x <= 180.0):
            return np.exp(-(180.0-x)/15)
        else:
            return 0
    
    hbond_score = fdist(x, bound, mu, sigma) * fangle(x, angle_cut)   
    
    return data
"""

#######################################################################################
################################### Initialisation ####################################
#######################################################################################

def fingerprint(protein, ligand, settings):
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
    
    bitstring = protein[2]
    
    fp1, bitstring = hbond_compute(protein[0], ligand[0], settings, bitstring)
    fp2, bitstring = saltbridge_compute(protein[0], ligand[0], settings, bitstring)
    fp3, bitstring = contact_compute(protein[0], ligand[0], settings, bitstring)
    fp4, bitstring = pistack_compute(protein[1], ligand[1], settings, bitstring)
    fp5, bitstring = pication_compute(protein[1], ligand[1], settings, bitstring)
    fp6, bitstring_ = pication_compute(ligand[1], protein[1], settings, bitstring)
    fp7, bitstring = Xbond_compute(protein[1], ligand[1], settings, bitstring)
    
    fp = fp1 + fp2 + fp3 + fp4 + fp5 + fp6 + fp7
    fp = pd.DataFrame.from_dict(fp)
    fp.insert(0, "index", int(protein[-1]), True)
    fp.insert(0, "ligand", ligand[-2], True)
    fp.insert(0, "protein", protein[-2], True)
    
    return fp, bitstring



settings = {}
settings['nonpolar'] = {6:1.7, 9:1.47, 17:1.75, 35:1.85, 53:1.98}
settings['hbond_dist_cut'] = 3.5
settings['hbond_angle_cut'] = 90.0

settings['hphob_dist_cut'] = 1.75
settings['hphob_dist_bound'] = 1.25

settings['contact_scale_cut'] = 1.75
settings['contact_scale_opt'] = 1.25

settings['saltbridge_resonance'] = True
settings['sbridge_dist_cut'] = 5.0

settings['pistack_dist_cut'] = 6.5
settings['pistack_angle_cut'] = 30.0

settings['pication_dist_cut'] = 6.5
settings['pication_angle_cut'] = 30.0

settings['Xbond_dist_cut'] = 5.0
settings['Xbond_angle_cut'] = 30.0

BASE_ANGLES = np.array((0, 180, 120, 109.5, 90), dtype=float)
    
fps = pd.DataFrame()


#######################################################################################
#################################### PATH DIRECTORY ###################################
#######################################################################################


trial = 'Trial05'
path = r'/home/yuyang/temp02'
coreset_pro = path + '/param/CoreSet_Receptor.txt'
coreset_lig = path + '/param/CoreSet_Ligand.txt'

with open(coreset_pro) as f:
    pro_list = [line.rstrip() for line in f]
    
with open(coreset_lig) as f:
    lig_list = [line.rstrip() for line in f]
    
    

#######################################################################################
################################### Main function #####################################
#######################################################################################


for i in range(len(pro_list)):
    pro_file = path + '/pocket/' + pro_list[i] + '_pocket.pdb'
    protein = MolFromPDBFile(pro_file)
    
    for j in range(4):
        k = 4 * i + j
        lig_files = glob.glob(path + '/vina/' + trial + '/output_native/' + 
                              lig_list[k] + '_out_model_' + '*.pdb')
        
        for l in range(len(lig_files)):
            ligand = MolFromPDBFile(lig_files[l])
            
            protein_addH = AddHs(protein, addCoords=True, addResidueInfo=True)
            ligand_addH = AddHs(ligand, addCoords=True)
            
            protein_ = Molecule(protein_addH, True, settings)
            ligand_ = Molecule(ligand_addH, False, settings)
            
            protein_oddt = next(oddt.toolkit.readfile('pdb', pro_file))
            protein_oddt.protein=True
            protein_oddt.addh()
            
            mols = list(oddt.toolkit.readfile('pdb', lig_files[l]))
            ligand_oddt = mols[0]
            ligand_oddt.addh()
            
            bitstring_pro = ['10000000']*len(protein_addH.GetAtoms())
            bitstring_lig = ['00000000']*len(ligand_addH.GetAtoms())
            
            
            """
            protein_ and ligand_ are RDKit Mol objects used to compute hydrogen
            bonds, salt bridges and hydrophobic contacts.
            
            protein_oddt and ligand_oddt are oddt.toolkit.Molecule objects used to
            compute pi stacking only.
            """
            
            protein_duo = [protein_, protein_oddt, bitstring_pro, pro_list[i], l]
            ligand_duo = [ligand_, ligand_oddt, bitstring_lig, lig_list[k], l]
            
            tmp, bitstring = fingerprint(protein_duo, ligand_duo, settings)
            
            fps = pd.concat([fps, tmp], sort=False)
            
            
            bitstring = [str(pro_list[i])] + bitstring 

            with open (trial + '_native_bitstring.txt', 'a') as fo:
                fo.write(','.join(str(z) for z in bitstring))
                fo.write('\n')


fps = fps[fps["protein_res"].str.contains("HOH") == False]
fps = fps[fps['dist'] <= 6.5] 

print(fps.head(10))

fps.to_csv(path + '/' + trial + '_native_descriptors.csv', index=False)



