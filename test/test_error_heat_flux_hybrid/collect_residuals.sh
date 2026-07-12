#!/bin/bash
# Usage: bash collect_residuals.sh SCHEME
# Outputs: residuals_SCHEME.dat with columns: fsize blsize residual

SCHEME=${1:-multi_point_iso}
OUTFILE="residuals_${SCHEME}.dat"
> $OUTFILE

for FSIZE in $(grep -v '^#\|^$' fsize.txt); do
  for BLSIZE in $(grep -v '^#\|^$' blsize.txt); do
    RESFILE="outputs/${SCHEME}/${FSIZE}_${BLSIZE}/residual.dat"
    if [ -f "$RESFILE" ]; then
      RESIDUAL=$(awk 'NR==1{r0=$3} END{print $3/r0}' $RESFILE)
    else
      RESIDUAL="NaN"
    fi
    echo "$FSIZE $BLSIZE $RESIDUAL" >> $OUTFILE
  done
  echo "" >> $OUTFILE
done

echo "Written: $OUTFILE"
