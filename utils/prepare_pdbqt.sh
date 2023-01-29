# This bash script aims to convert all receptor and ligand files in .pdb to .pdbqt format.

cd /home/yuyang/temp02/param

for (( i=1; i<=57; i++ )) 
do    
    
    cd /home/yuyang/temp02/param
    PRO_HANDLE=$(cat CoreSet_Receptor.txt | head -$i | tail -1)
    
    cd /home/yuyang/temp02/protein
    python prepare_receptor4.py -l ${PRO_HANDLE}_protein.pdb -v -o ${PRO_HANDLE}.pdbqt -A hydrogens
    
    mv *.pdbqt ../protein_pdbqt
    
done

for (( j=1; j<=228; j++ )) 
do    
    
    cd /home/yuyang/temp02/param
    LIG_HANDLE=$(cat CoreSet_Ligand.txt | head -$j | tail -1)
    
    cd /home/yuyang/temp02/ligand
    python prepare_ligand4.py -l ${LIG_HANDLE}_ligand.pdb -v -o ${LIG_HANDLE}.pdbqt -A hydrogens
    
    mv *.pdbqt ../ligand_pdbqt
    
done

exit;

