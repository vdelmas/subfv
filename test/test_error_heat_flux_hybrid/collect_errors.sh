#!/bin/bash
# Usage: bash collect_errors.sh SCHEME
# Outputs: errors_SCHEME.dat with columns: fsize blsize err_st err_cp

SCHEME=${1:-multi_point_iso}
OUTFILE="errors_${SCHEME}.dat"
> $OUTFILE

for FSIZE in $(grep -v '^#\|^$' fsize.txt); do
  for BLSIZE in $(grep -v '^#\|^$' blsize.txt); do
    ERRFILE="outputs/${SCHEME}/${FSIZE}_${BLSIZE}/error.txt"
    if [ -f "$ERRFILE" ]; then
      ERROR=$(cat $ERRFILE)
    else
      ERROR="NaN NaN"
    fi
    echo "$FSIZE $BLSIZE $ERROR" >> $OUTFILE
  done
  echo "" >> $OUTFILE
done

echo "Written: $OUTFILE"
