# The aim of this script is to read Trial**_rmsd_bool.txt
# If the bool value is 2, then it is a native case. The respective ligand model will be moved to output_native
# If the bool value is 3, then it is a nonnative case. The respective ligand will be moved to output_refernece


XPATH=/home/yuyang/temp02

LINE=0

for ((i=1; i<=228; i++))
do
    cd $XPATH/param
    LIG_HANDLE=$(cat CoreSet_Ligand_Sorted.txt | head -$i | tail -1)
    
    cd $XPATH/vina/Trial05/output_split
    NUM_FILES=$(grep -l "$LIG_HANDLE" * | wc -l)
    
    for ((j=1; j<=${NUM_FILES}; j++))
    do
        cd $XPATH/vina
        LINE=$(( ( $i - 1 ) * 20 ))
        LINE=$(( $LINE + $j ))
        BOOL=$(cat Trial05_rmsd_bool.txt | head -$LINE | tail -1)
        VALUE=$(echo ${BOOL:0:1})
        
        if (( $VALUE == 2 )); then
            cd $XPATH/vina/Trial05/output_split
            cp ${LIG_HANDLE}_out_model_$j.pdb ../output_native
        elif (( $VALUE == 3 )); then
            cd $XPATH/vina/Trial05/output_split
            cp ${LIG_HANDLE}_out_model_$j.pdb ../output_reference
        else
            echo "Rejected:" ${LIG_HANDLE} $j
        fi
        
    done   
done



cd $XPATH/vina/Trial05/output_split
NUM_FILES_ORI=$(ls . | wc -l)
echo $NUM_FILES_ORI


cd $XPATH/vina/Trial05/output_native
NUM_FILES_NVE=$(ls . | wc -l)
echo $NUM_FILES_NVE


cd $XPATH/vina/Trial05/output_reference
NUM_FILES_REF=$(ls . | wc -l)
echo $NUM_FILES_REF

exit;
