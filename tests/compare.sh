#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

#actual output metrics is $1, expected output metrics is $2
MATCHED=$(join --nocheck-order $1 $2 | wc -l)
NOMATCH1=$(join --nocheck-order -v 1 $1 $2 | wc -l)
NOMATCH2=$(join --nocheck-order -v 2 $1 $2 | wc -l)
ALL=$(( $NOMATCH1 + $NOMATCH2 + $MATCHED))
MATCH_RATIO=$(( $MATCHED * 100 / $ALL ))
if [[ $MATCH_RATIO>=90 ]];then
  exit 0
else
  echo "Only $MATCH_RATIO% variants match, FAILED"
fi

exit 1
~
