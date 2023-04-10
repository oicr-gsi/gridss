#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#actual output metrics is $1, expected output metrics is $2
MATCHED=$(cat $1 $2 | sort -V | uniq -d | wc -l)
NOMATCH=$(cat $1 $2 | sort -V | uniq -u | wc -l)
ALL=$(( $NOMATCH + $MATCHED))
MATCH_RATIO=$(( $MATCHED * 100 / $ALL ))
echo "We have $MATCH_RATIO% variants matched"
if [ $MATCH_RATIO -gt 90 ] ; then
  exit 0
else
  echo "Only $MATCH_RATIO% variants match, FAILED"
fi
~
