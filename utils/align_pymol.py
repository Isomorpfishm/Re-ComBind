import __main__
import glob
import os
import pymol
import time



def align(mobile, reference, filename, resname):
    
    """Do a sequence alignment between two pdbs and write the output.
    Parameters
    ----------
    mobile : A PDB file. Sequence will be first in the alignment.
    reference : A PDB file. Sequence will be second in alignment.
    filename : Output name of alignment file.
    """
    __main__.pymol_argv = ['pymol', '-qc']
    
    pymol.finish_launching()
    pymol.cmd.load(mobile)
    pymol.cmd.load(reference)
    obj = list(pymol.cmd.get_object_list('all'))
    pymol.cmd.align(obj[0], obj[1], object='pdb', transform=1)
    pymol.cmd.select('sele', 'resname ' + resname)
    pymol.cmd.save(filename, selection = 'sele', format = 'pdb')
    
    com = pymol.cmd.centerofmass('sele')
    
    pymol.cmd.reinitialize()
    time.sleep(1)
    
    return com
    
    
# Specify file directory
path = r'/home/yuyang/temp'
coreset_pro = r'/home/yuyang/temp/param/CoreSet_Receptor.txt'
coreset_lig = r'/home/yuyang/temp/param/CoreSet_Ligand.txt'
lig_resname = r'/home/yuyang/temp/param/Ligand_resname.txt'

with open(coreset_pro) as f:
    ref_list = [line.rstrip() for line in f]
    
with open(coreset_lig) as f:
    mobile_list = [line.rstrip() for line in f]
    
with open(lig_resname) as f:
    resname_list = [line.rstrip() for line in f]

    
# Initialise center of mass (com) list
com = []

for i in range(len(ref_list)):
    for j in range(4):
            k = 4 * i + j
            
            ref_pdb = path + '/rigid/' + ref_list[i] + '.pdb'
            mobile_pdb = path + '/rigid/' + mobile_list[k] + '.pdb'
            
            output = path + '/aligned/' + mobile_list[k] + '_aligned.pdb'
            com_tmp = align(mobile = mobile_pdb, reference = ref_pdb, filename = output, resname = resname_list[k])
            
            com.append(com_tmp)
            

# Save the center of mass to a .txt file
with open("centerofmass.txt", "a") as f:
    for item in com:
        f.write(str(item) + "\n")
