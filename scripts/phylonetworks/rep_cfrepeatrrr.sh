#!/bin/bash
# shell script to run SNPs2CF on multiple replicates of snp datasets.
# 
# author: thomas markus huber
# contact: thomas.huber@evobio.eu / thomasmarkus.huber.3696@student.uu.se
# input needed: - directory with rep-folders with .phy and imap.txt files
# 		  (made with fastapasta_snphywatzko.py)

echo "~~~ Start of SNPs2CF runs for replicates within $1 ... ~~~"
PHY="$1/pnw_data.phy"
IMAP="$1/pnw_imap.txt"
echo "  - run for $PHY and $IMAP in $1 started ..."
R --vanilla --slave --args $PHY $IMAP < ./scripts/phylonetworks/ms_concordiafactoria.r
echo "  ... run for $PHY and $IMAP in $1 done -"


# R --vanilla --slave --args INPUT < ./scripts/phylonetworks/ms_concordiafactoria.r
