#!/bin/bash
# shell script for the HyDe phylogenetic network analysis
# author: thomas markus huber
# contact: thomas.huber@evobio.eu / thomasmarkus.huber.3696@student.uu.se
# input needed: - phylip formatted sequence file
#		- map txt file linking sequence/indiv names with species/OTUs (tab-delim)
#		- triples file listing the individuals to compare (parent1, hybrid, parent2)
#			can be output of run_hyde.py

echo "~~~~~~~~~~ This is the trial of HyDe for phylogenetic networks ~~~~~~~~~~"
BIGDIR=$1

# INPUT="`sed '8q;d' $INPUT_TXT`"
# MAP="${BIGDIR}/2020_04_vitis_USA-map.txt"
MAP="${BIGDIR}/vitisasia-map.txt"
# OUT="vsylveast"
OUT="vsylvwest"
# TRIP="${BIGDIR}/wdwreps-triples.txt"
TRIP="${BIGDIR}/vitisasia-triples.txt"

# OUTPUT_FILE="`sed '11q;d' $INPUT_TXT`"

NUM_IND=`awk '{ print $1 }' $MAP | sort | uniq | awk 'END {print NR}'`
NUM_TAXA=`awk '{ print $2 }' $MAP | sort | uniq | awk 'END {print NR}'`


for REP in `ls "${BIGDIR}/data"`
do
	echo $REP
    INPUT="${BIGDIR}/data/${REP}"
    OUTPUT_FILE="${BIGDIR}/out/${REP}"
	echo "~ input-file: $INPUT ~ map-file: $MAP ~ triples-file: $TRIP ~ $NUM_IND individuals in $NUM_TAXA with $NUM_SITES ~ outgroup: $OUT ~"
	NUM_SITES=$((`awk '{ print $2 }' $INPUT | head -1 | wc -c`-1))

	individual_hyde.py -i $INPUT -m $MAP -tr $TRIP -o $OUT -n $NUM_IND -t $NUM_TAXA -s $NUM_SITES --prefix $OUTPUT_FILE
        # individual_hyde.py -i $INPUT -m $MAP -tr $TRIP -o $OUT -n $NUM_IND -t $NUM_TAXA -s $NUM_SITES

	# echo "~~~~~~~~~~ individual_hyde is done and saved as/in ${OUTPUT_FILE}-ind.txt ~~~~~~~~~~"
	# Filter (not needed):
	# head -1 ${OUTPUT_FILE}-ind.txt > ${OUTPUT_FILE}-ind-filtered.txt
	# grep . ${OUTPUT_FILE}-ind.txt | awk '$5<0.05' >> ${OUTPUT_FILE}-ind-filtered.txt

	# echo "~~~~~~~~~~ an additional bonferroni correction of the p_values will be performed with and r_bonferroni.r ~~~~~~~~~~"
	R --vanilla --slave --args ${OUTPUT_FILE}-ind.txt < ./scripts/hyde/r_bonferroni.r
	# Filter:
	# head -1 ${OUTPUT_FILE}-ind-bonf.txt > ${OUTPUT_FILE}-ind-bonfiltered.txt
	# grep . ${OUTPUT_FILE}-ind-bonf.txt | awk '$5<1.0' >> ${OUTPUT_FILE}-ind-bonfiltered.txt
	# echo "~~~~~~~~~~ done. the whole set was save as ${OUTPUT_FILE}-ind-bonf.txt and filtered by significance as ${OUTPUT_FILE}-ind-bonf.txt ~~~~~~~~~~"

done
