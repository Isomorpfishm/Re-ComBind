import oddt
from oddt import interactions

import numpy as np

import glob



# Specify path directory and protein-ligand sources
trial = 'Trial01'
path = r'/home/yuyang/temp02'
coreset_pro = path + '/param/CoreSet_Receptor.txt'
coreset_lig = path + '/param/CoreSet_Ligand.txt'

with open(coreset_pro) as f:
    rep_list = [line.rstrip() for line in f]
    
with open(coreset_lig) as f:
    lig_list = [line.rstrip() for line in f]
    
    
    
# Finding interactions using the ODDT library
class find_interactions:
    def __init__(self, protein, ligand):
        self.protein = protein
        self.ligand = ligand
        
    def hbonds(self):
        acceptor_atoms, donor_atoms, strict = interactions.hbond_acceptor_donor(self.protein, self.ligand)
        
        acc_id, don_id = acceptor_atoms['id'].tolist(), donor_atoms['id'].tolist()
        acc_coords, don_coords = acceptor_atoms['coords'].tolist(), donor_atoms['coords'].tolist()
        acc_resname, don_resname = acceptor_atoms['resname'].tolist(), donor_atoms['resname'].tolist()
        
        protein_vec, ligand_vec = [], []
        
        for p in range(len(protein.atoms)):
            protein_vec.append(0)
                
        for q in range(len(ligand.atoms)):
            ligand_vec.append(0)
            
        for idx in acc_id:
            protein_vec[idx-1] = 1
            
        for idy in don_id:
            ligand_vec[idy-1] = 2
        
        hbond_dict = {}
        hbond_dict["acceptor_atoms_resname"] = acc_resname
        hbond_dict["donor_atoms_resname"] = don_resname
        hbond_dict["protein_vector"] = protein_vec
        hbond_dict["ligand_vector"] = ligand_vec
        
        # HBONDS: Check an atom belongs to protein or ligand
        acceptor_source = check_aa(hbond_dict["acceptor_atoms_resname"])
        donor_source = check_aa(hbond_dict["donor_atoms_resname"])
        
        acceptor_ele, acceptor_coords, donor_ele, donor_coords = check_hydrogen(acceptor_atoms, donor_atoms, self.protein, self.ligand, donor_source)
        
        hbond_dict["D-H_element"] = donor_ele
        hbond_dict["D-H_coords"] = donor_coords
        hbond_dict["A_element"] = acceptor_ele
        hbond_dict["A_coords"] = acceptor_coords
        
        return hbond_dict
        
    def hphobs(self):
        pro, lig = interactions.hydrophobic_contacts(self.protein, self.ligand)
        
        pro_id, lig_id = pro['id'].tolist(), lig['id'].tolist()
        pro_r, lig_r = pro['radius'].tolist(), lig['radius'].tolist()
        
        for p in range(len(pro_id)):
            pro_id[p] += 1
            
        for q in range(len(lig_id)):
            lig_id[q] += 1
            
        pro_res = pro_id + pro_r
        lig_res = lig_id + lig_r
         
        hphob_dict = {}
        hphob_dict["protein_id_radius"] = pro_res
        hphob_dict["ligand_id_radius"] = lig_res
                    
        return hphob_dict
        
    def sbridges(self):
        plus_atoms, minus_atoms = interactions.salt_bridge_plus_minus(self.protein, self.ligand)
        
        plus_atoms_id, minus_atoms_id = plus_atoms['id'].tolist(), minus_atoms['id'].tolist()
        plus_atoms_r, minus_atoms_r = plus_atoms['radius'].tolist(), minus_atoms['radius'].tolist()
        
        for p in range(len(plus_atoms_id)):
            plus_atoms_id[p] += 1
            
        for q in range(len(minus_atoms_id)):
            minus_atoms_id[q] += 1
            
        plus_atoms_res = plus_atoms_id + plus_atoms_r
        minus_atoms_res = minus_atoms_id + minus_atoms_r
         
        sbridge_dict = {}
        sbridge_dict["plus_atoms_id_radius"] = plus_atoms_res
        sbridge_dict["minus_atoms_id_radius"] = minus_atoms_res
                    
        return sbridge_dict
        
    def pistack(self):
        ring1, ring2, strict_para, strict_perp = interactions.pi_stacking(self.protein, self.ligand)
        
        ring1_ctrid, ring2_ctrid = ring1['centroid'].tolist(), ring2['centroid'].tolist()
        ring1_vec, ring2_vec = ring1['vector'].tolist(), ring2['vector'].tolist()
        ring1_resnum, ring1_resname = ring1['resnum'].tolist(), ring2['resname'].tolist()
        
        rings_ctr, rings_vec, ring1_res = [], [], []
        
        rings_ctr = ring1_ctrid + ring2_ctrid
        rings_vec = ring1_vec + ring2_vec
        
        for p in range(len(ring1_resnum)):
            ring1_res.append(ring1_resname[p])
            ring1_res.append(ring1_resnum[p]) 
        
        pistack_dict = {}
        pistack_dict["rings_centroid"] = rings_ctr
        pistack_dict["rings_vector"] = rings_vec
        pistack_dict["ring1_residue"] = ring1_res
        
        
        return pistack_dict

    

# Check if an atom belongs to protein or ligand based on resname
# If resname == '' then it comes from the ligand
# Else it comes from the protein
def check_aa(input_list):
    source = []
    AA_LIST = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN',
               'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
               'MET', 'PHE', 'PRO', 'PYL', 'SER', 'SEC',
               'THR', 'TRP', 'TYR', 'VAL', 'ASX', 'GLX']
    
    for residue in input_list:
        if residue in AA_LIST:
            source.append('PROTEIN')
        else:
            source.append('LIGAND')
            
    return source


# Find the donor-hydrogen elements and their coordinates
def check_hydrogen(acceptor_list, donor_list, protein, ligand, source):
    don_id = donor_list['id'].tolist()
    don_neighbors_id = donor_list['neighbors_id'].tolist()
    don_coords = donor_list['coords'].tolist()
    don_ele = donor_list['atomicnum'].tolist()
    
    acc_id = acceptor_list['id'].tolist()
    acc_neighbors_id = acceptor_list['neighbors_id'].tolist()
    acc_coords = acceptor_list['coords'].tolist()
    acc_ele = acceptor_list['atomicnum'].tolist()
    
    acceptor_ele, acceptor_coords, donor_ele, donor_coords = [], [], [], []
    
    assert len(donor_list) == len(source)
    
    for i in range(len(source)):
        cnt = 0
        for j in range(6):
            if source[i] == 'PROTEIN':
                if don_neighbors_id[i][j] != 0:
                    ngh_ele = protein.atoms[don_neighbors_id[i][j]].atomicnum
                    
                    if ngh_ele == 1:
                        cnt += 1
                        
                        donor_ele.append(don_ele[i])
                        donor_coords.append(don_coords[i])
                        donor_ele.append(1)
                        h_coords = protein.atoms[don_neighbors_id[i][j]].coords
                        donor_coords.append(h_coords)
                        
                        acceptor_ele.append(acc_ele[i])
                        acceptor_coords.append(acc_coords[i])
                    
                    else:
                        pass
                
                else:
                    pass
                
            elif source[i] == 'LIGAND':
                if don_neighbors_id[i][j] != 0:
                    ngh_ele = ligand.atoms[don_neighbors_id[i][j]].atomicnum
                    
                    if ngh_ele == 1:
                        cnt += 1
                        
                        donor_ele.append(don_ele[i])
                        donor_coords.append(don_coords[i])
                        donor_ele.append(1)
                        h_coords = ligand.atoms[don_neighbors_id[i][j]].coords
                        donor_coords.append(h_coords)
                        
                        acceptor_ele.append(acc_ele[i])
                        acceptor_coords.append(acc_coords[i])
                    
                    else:
                        pass
                
                else:
                    pass
                    
            else:
                pass
                
    
    return acceptor_ele, acceptor_coords, donor_ele, donor_coords
        


# Appending protein and ligand information to the output
def mod_list(result_list, num, rep, lig):
    result_list.insert(0, int(num))
    result_list.insert(0, str(lig))
    result_list.insert(0, str(rep))
    
    return result_list
    

# Export results
def export_txt(filename, result_list):
    filehandle = open(filename, 'a')
    print(*result_list, sep = ",", file = filehandle)
    filehandle.close()

    


# Main function
for i in range(len(rep_list)):
    # Index i controls the receptors (there are 57 of them)
    
    rep_files = path + '/pocket/' + rep_list[i] + '_pocket.pdb'
    protein = next(oddt.toolkit.readfile('pdb', rep_files))
    protein.protein = True
    protein.addh()
    
    for j in range(4):
        # Index k controls the ligands (there are 228 of them)
        # Each receptor has 4 ligands
        
        k = 4 * i + j

        lig_files = glob.glob(path + '/vina/' + trial + '/output_filter/' + lig_list[k] + '_out_model_' + '*.pdb')
        
        for l in range(len(lig_files)):
            # Index l controls the number of models for each ligand
            # E.g. there are 16 mols for ligand 1eby
            # So index l goes up to 16
            
            mols = list(oddt.toolkit.readfile('pdb', lig_files[l]))
            ligand = mols[0]
            ligand.addh()
            
            all_interactions = find_interactions(protein, ligand)
            
            # Finding hydrogen bonds 
            hbond_result = all_interactions.hbonds()
            
            # Finding hydrophobic contacts
            hphob_result = all_interactions.hphobs()
            
            # Finding salt bridges
            sbridge_result = all_interactions.sbridges()
            
            # Finding pi_stacking
            pistack_result = all_interactions.pistack()
            
            
            
            
            #acceptor = acceptor_source + hbond_result["acceptor_atoms_coords"]
            #donor = donor_source + hbond_result["donor_atoms_coords"]
            
            
            # Appending protein and ligand information to the output
            hbond_result["protein_vector"] = mod_list(result_list = hbond_result["protein_vector"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            hbond_result["ligand_vector"] = mod_list(result_list = hbond_result["ligand_vector"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            hbond_result["D-H_element"] = mod_list(result_list = hbond_result["D-H_element"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            hbond_result["D-H_coords"] = mod_list(result_list = hbond_result["D-H_coords"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            hbond_result["A_element"] = mod_list(result_list = hbond_result["A_element"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            hbond_result["A_coords"] = mod_list(result_list = hbond_result["A_coords"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            
            
            #hphob_result["protein_id_radius"] = mod_list(result_list = hphob_result["protein_id_radius"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            #hphob_result["ligand_id_radius"] = mod_list(result_list = hphob_result["ligand_id_radius"], num = l+1, lig = lig_list[k], rep = rep_list[i])

            #sbridge_result["plus_atoms_id_radius"] = mod_list(result_list = sbridge_result["plus_atoms_id_radius"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            #sbridge_result["minus_atoms_id_radius"] = mod_list(result_list = sbridge_result["minus_atoms_id_radius"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            
            #pistack_result["rings_centroid"] = mod_list(result_list = pistack_result["rings_centroid"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            #pistack_result["rings_vector"] = mod_list(result_list = pistack_result["rings_vector"], num = l+1, lig = lig_list[k], rep = rep_list[i])
            #pistack_result["ring1_residue"] = mod_list(result_list = pistack_result["ring1_residue"], num = l+1, lig = lig_list[k], rep = rep_list[i])

            
            # Export results
            export_txt(filename = path + '/vina/' + trial + '_hbond_pro_vector.txt', result_list = hbond_result["protein_vector"])
            export_txt(filename = path + '/vina/' + trial + '_hbond_lig_vector.txt', result_list = hbond_result["ligand_vector"])
            export_txt(filename = path + '/vina/' + trial + '_hbond_D-H_element.txt', result_list = hbond_result["D-H_element"])
            export_txt(filename = path + '/vina/' + trial + '_hbond_D-H_coords.txt', result_list = hbond_result["D-H_coords"])
            export_txt(filename = path + '/vina/' + trial + '_hbond_A_elements.txt', result_list = hbond_result["A_element"])
            export_txt(filename = path + '/vina/' + trial + '_hbond_A_coords.txt', result_list = hbond_result["A_coords"])                 
            #export_txt(filename = path + '/vina/' + trial + '_hphob_pro_id_radius.txt', result_list = hphob_result["protein_id_radius"])
            #export_txt(filename = path + '/vina/' + trial + '_hphob_lig_id_radius.txt', result_list = hphob_result["ligand_id_radius"])
            #export_txt(filename = path + '/vina/' + trial + '_sbridge_plus_id_radius.txt', result_list = sbridge_result["plus_atoms_id_radius"])
            #export_txt(filename = path + '/vina/' + trial + '_sbridge_minus_id_radius.txt', result_list = sbridge_result["minus_atoms_id_radius"])
            #export_txt(filename = path + '/vina/' + trial + '_pistack_rings_centroid.txt', result_list = pistack_result["rings_centroid"])
            #export_txt(filename = path + '/vina/' + trial + '_pistack_rings_vector.txt', result_list = pistack_result["rings_vector"])
            #export_txt(filename = path + '/vina/' + trial + '_pistack_ring1_residue.txt', result_list = pistack_result["ring1_residue"])   
                       
