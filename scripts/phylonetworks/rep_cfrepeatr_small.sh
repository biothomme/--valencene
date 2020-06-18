#!/bin/bash
# shell script to run SNPs2CF on multiple replicates of snp datasets.
# split up run 16-20
# author: thomas markus huber
# contact: thomas.huber@evobio.eu / thomasmarkus.huber.3696@student.uu.se
# input needed: - directory with rep-folders with .phy and imap.txt files
# 		  (made with fastapasta_snphywatzko.py)

echo "~~~ Start of SNPs2CF runs for replicates within $1 ... ~~~"
BIGDIR=`ls $1`
for DIR in $BIGDIR
do
    REPNUM=`echo $DIR | grep -E -o '[1-9]+[0-9]?+'`
    if [ $(($REPNUM)) -gt 15 ] && [ $(($REPNUM)) -lt 21 ]
    then
        PHY="$1/$DIR/pnw_data.phy"
        IMAP="$1/$DIR/pnw_imap.txt"
        echo "  - run for $PHY and $IMAP started ..."
        R --vanilla --slave --args $PHY $IMAP < ./scripts/phylonetworks/ms_concordiafactoria.r
        echo "  ... run for $PHY and $IMAP done -"
    fi
done


