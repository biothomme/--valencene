#!/bin/bash
#
# this script makes replicates of the dataset with different amounts of ibs clusters
# input: directory, to save reps

IBSDIR="data/real_data/vitis/additional_information/ibs_cluster"
SMRY="data/real_data/vitis/2020_04_vitis_USA_2_summary.csv"
PERMUCLUST="scripts/data_analysis/vino_permuclust.R"
OUTPUT=$1


for i in {2..16}
do
	CS=$(($i ** (2) / 4))
        # FOLD="${OUTPUT}/vs_${CS}/"
	# mkdir $FOLD
	for j in {2..16}
	do
		CV=$(($j ** (2) / 4))
		FOLD="${OUTPUT}/vs_${CS}_vv_${CV}/"
		echo $FOLD
		mkdir $FOLD
		R --vanilla --slave --args $IBSDIR $SMRY $CV $CS $FOLD < $PERMUCLUST 
	done
done

