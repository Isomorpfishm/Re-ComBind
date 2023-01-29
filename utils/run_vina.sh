# This bash script is to run QuickVina 2 and organise the output log (.txt) and ligand (.pdbqt) files.


QVINA=/home/ygmu/softwares/qvina/bin/qvina02
PTH=/home/yuyang/temp02

for ((i=1; i<=57; i++))
do
	cd $PTH/param

	PRO_HANDLE=$(cat CoreSet_Receptor.txt | head -$i | tail -1)
	PRO_NAME=$PTH/protein_pdbqt/${PRO_HANDLE}.pdbqt

	echo "Protein file name is $PRO_NAME"

    cd $PTH/vina/Trial05
	mkdir Protein_${PRO_HANDLE}

	for ((j=1; j<=4; j++))
	do
		count=$((4 * ($i - 1) + $j))

		cd $PTH/param

		LIG_HANDLE=$(cat CoreSet_Ligand.txt | head -$count | tail -1)
		LIG_NAME=$PTH/ligand_pdbqt/${LIG_HANDLE}.pdbqt

		CMX=$(cat com_x.txt | head -$count | tail -1)
		CMY=$(cat com_y.txt | head -$count | tail -1)
		CMZ=$(cat com_z.txt | head -$count | tail -1)

		echo "Coordinates of center of mass:" "$CMX" "," "$CMY" "," "$CMZ"

		cd $PTH/vina/Trial05/Protein_${PRO_HANDLE}

		$QVINA --receptor $PRO_NAME --ligand $LIG_NAME --center_x $CMX --center_y $CMY --center_z $CMZ --size_x 15 --size_y 15 --size_z 15 --cpu 24 --num_modes 20 > $LIG_HANDLE.txt

		cd $PTH/ligand_pdbqt
		mv ${LIG_HANDLE}_out.pdbqt $PTH/vina/Trial05/Protein_${PRO_HANDLE}

		echo "Processed $i out of 57 proteins WITH $j out 4 of ligands."
	done
done

exit; 
