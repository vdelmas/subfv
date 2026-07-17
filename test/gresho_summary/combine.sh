#!/bin/bash
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
SCHEMES_FILE="$ROOT/../schemes.txt"

PDFS=()
while IFS= read -r SCHEME; do
  PDF="$ROOT/outputs/$SCHEME/summary.pdf"
  if [ -f "$PDF" ]; then
    PDFS+=("$PDF")
  else
    echo "Warning: $PDF not found, skipping" >&2
  fi
done < "$SCHEMES_FILE"

pdfunite "${PDFS[@]}" "$ROOT/outputs/gresho_summary_raw.pdf"

TMPIN=$(mktemp /tmp/gresho_raw_XXXXXX.pdf)
cp "$ROOT/outputs/gresho_summary_raw.pdf" "$TMPIN"

TMPOUT=$(mktemp /tmp/gresho_summary_XXXXXX.pdf)
gs -q -dBATCH -dNOPAUSE -sDEVICE=pdfwrite \
   -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer \
   -sOutputFile="$TMPOUT" "$TMPIN"
mv "$TMPOUT" "$ROOT/outputs/gresho_summary.pdf"

TMPOUT=$(mktemp /tmp/gresho_summary_XXXXXX.pdf)
gs -q -dBATCH -dNOPAUSE -sDEVICE=pdfwrite \
   -dCompatibilityLevel=1.4 -dPDFSETTINGS=/ebook \
   -sOutputFile="$TMPOUT" "$TMPIN"
mv "$TMPOUT" "$ROOT/outputs/gresho_summary_light.pdf"

rm "$TMPIN" "$ROOT/outputs/gresho_summary_raw.pdf"

date > "$ROOT/outputs/combine.timestamp"
