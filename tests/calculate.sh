#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#enter the workflow's final output directory ($1)
cd $1

#find all vcf files, return their PASS call ids to std out
for v in *.vcf;do cat $v | grep PASS | awk '{OFS="_";print $1,$2,$4,$5}';done 
