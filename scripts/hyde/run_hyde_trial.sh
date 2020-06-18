#!/bin/bash
# shell script for the HyDe phylogenetic network analysis
# author: thomas markus huber
# contact: thomas.huber@evobio.eu / thomasmarkus.huber.3696@student.uu.se
# input needed: - phylip formatted sequence file
#		- map txt file linking sequence/indiv names with species/OTUs (tab-delim)


echo "~~~~~~~~~~ This is the trial of HyDe for phylogenetic networks ~~~~~~~~~~"

INPUT_TXT=$1
INPUT="`sed '7q;d' $INPUT_TXT`"
MAP="`sed '8q;d' $INPUT_TXT`"
OUT="`sed '9q;d' $INPUT_TXT`"
OUTPUT_FILE="`sed '10q;d' $INPUT_TXT`"
NUM_IND=`awk '{ print $1 }' "$MAP" | sort | uniq | awk 'END {print NR}'`
NUM_TAXA=`awk '{ print $2 }' "$MAP" | sort | uniq | awk 'END {print NR}'`
NUM_SITES=$((`awk '{ print $2 }' "$INPUT" | head -1 | wc -c`-1))

echo "~ input-file: $INPUT ~ map-file: $MAP ~ $NUM_IND individuals in $NUM_TAXA taxa with $NUM_SITES loci ~ outgroup: $OUT ~"

run_hyde.py -i $INPUT -m $MAP -o $OUT -n $NUM_IND -t $NUM_TAXA -s $NUM_SITES --prefix $OUTPUT_FILE

echo "~~~~~~~~~~ run_hyde is done and saved as/in ${OUTPUT_FILE}-out.txt ~~~~~~~~~~"
# Filter:
head -1 ${OUTPUT_FILE}-out.txt > ${OUTPUT_FILE}-out-filtered.txt
grep . ${OUTPUT_FILE}-out.txt | awk '$5<0.05' >> ${OUTPUT_FILE}-out-filtered.txt


echo "~~~~~~~~~~ an additional bonferroni correction of the p_values will be performed with and r_bonferroni.r ~~~~~~~~~~"
R --vanilla --slave --args ${OUTPUT_FILE}-out.txt < ./scripts/hyde/r_bonferroni.r
# Filter:
head -1 ${OUTPUT_FILE}-out-bonf.txt > ${OUTPUT_FILE}-out-bonfiltered.txt
grep . ${OUTPUT_FILE}-out-bonf.txt | awk '$5<1.0' >> ${OUTPUT_FILE}-out-bonfiltered.txt
echo "~~~~~~~~~~ done. the whole set was save as ${OUTPUT_FILE}-out-bonf.txt and filtered by significance as ${OUTPUT_FILE}-out-bonf.txt ~~~~~~~~~~"



