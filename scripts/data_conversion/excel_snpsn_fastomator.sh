#!/bin/bash
# shell script to transform excel snp sheets (like roberto bacilieris grapes)
# to a smooth fasta file
# author: thomas markus huber
# contact: thomas.huber@evobio.eu / thomasmarkus.huber.3696@student.uu.se
# input needed: - excel sheet with snp data labelled by populations 
#		  and individuals
# 
# the script randomly picks one of the two diploid alleles and thus
# pseudo-haploidizes the sequences.
# therefore $RANDOM is used which produces a random number between 0 and 32767.
# the lower 50 % are from 0 to 16383, the upper from 16384 to 32767.
# in addition: NAs are excluded with if condition, they were dangerous - could lead to base A!

echo "~~~ Start of transformation of $1 into fasta format ... ~~~"

# this needs to be done for robertos files:
dos2unix $1

LENGTH=`head -1 $1 | awk '{ print NF }'`
SEQNUM=$((`awk 'END {print NR}' $1`-3))

touch $1.fasta

for i in $(seq 2 1 $LENGTH)
do
	POP=`sed '3q;d' $1 | tr -d ' ' | cut -f$i`
	ADN=`cut -f$i $1 | sed '1q;d'`

        HEADER=`echo ">${POP}|${ADN}" | tr -d "\r"`
	BASES1=`cut -f$i $1 | tail -n $SEQNUM | cut -c1-1 | tr -d '\n'`
	BASES2=`cut -f$i $1 | tail -n $SEQNUM | cut -c2-2 | tr -d '\n'`
	BASES=""
	for j in $(seq 2 1 $(($SEQNUM+1)))
	do
                # the following commented if comparison can be used to avoid NAs, if they are present
		# if [[ "${BASES1:$j:1}" == "N" ]]
		# then
		# 	BASES="${BASES}${BASES1:$j:1}"
		# else
        if (( $RANDOM > 16383 ))
		then
			BASES="${BASES}${BASES1:$j:1}"
		else
            BASES="${BASES}${BASES2:$j:1}"
		fi
		# fi
	done
	echo "${HEADER}" >> $1.fasta
    echo "${BASES}" >> $1.fasta
done

echo "~~~ File was transformed into fasta format and saved as $1.fasta. ~~~"


