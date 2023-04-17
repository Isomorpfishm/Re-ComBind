#!/bin/bash


while getopts 'p:l:' OPTION; do
    case "$OPTION" in
        p)
            PRO_FOL="$OPTARG"
            ;;
        l)
            LIG_FOL="$OPTARG"
            ;;
        ?)
            echo "script usage: $(basename \$0) [-p folderOfProtein][-l folderOfLigand] [-n numOfDocking]" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
          ;;
    esac
done
shift "$(($OPTIND -1))"

# Define the folder name
folder_name="descriptors"

# Check if the folder exists
if [ ! -d "$folder_name" ]; then
  # If the folder does not exist, create it
  mkdir $folder_name
fi


python ifp.py -p $PRO_FOL -l $LIG_FOL
