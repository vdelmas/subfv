#!/bin/bash
# gen_sum_pdf.sh — Generate a summary PDF for one scheme across multiple test cases.
# Usage: ./gen_sum_pdf.sh <scheme> [test1 test2 ...]
# If no tests are listed, all test directories are scanned.
# Output: <scheme>_summary.pdf in the current directory.

set -euo pipefail

SCHEME=${1:?Usage: $0 <scheme> [test1 test2 ...]}
shift

TEST_DIR=$(cd "$(dirname "$0")"; pwd)
OUT_PDF="$(pwd)/${SCHEME}_summary.pdf"
WORK=$(mktemp -d)
trap "rm -rf $WORK" EXIT
TEX="$WORK/summary.tex"

# ---------- helpers ----------------------------------------------------------

ltx() { echo "$1" | sed 's/_/\\_/g; s/&/\\&/g'; }

minipage() {
    local w=$1 local=$2 caption=$3
    cat >> "$TEX" <<FIG
  \begin{minipage}[t]{${w}\linewidth}\centering
    \includegraphics[width=\linewidth]{${local}}\\\\[2pt]
    {\small ${caption}}
  \end{minipage}\hfill
FIG
}

# ---------- test list --------------------------------------------------------

if [ $# -eq 0 ]; then
    _all=$(find "$TEST_DIR" -mindepth 1 -maxdepth 1 -type d | xargs -n1 basename | sort)
    mapfile -t TESTS < <(
        echo "$_all" | grep '^toro'        ;
        echo "$_all" | grep '^test_shear'  ;
        echo "$_all" | grep '^gresho'      ;
        echo "$_all" | grep '^sedov'       ;
        echo "$_all" | grep '_euler$'      ;
        echo "$_all" | grep '_ns$'         ;
        echo "$_all" | grep -vE '^toro|^test_shear|^gresho|^sedov|_euler$|_ns$'
    )
else
    TESTS=("$@")
fi

# ---------- LaTeX header -----------------------------------------------------

cat > "$TEX" <<HEADER
\documentclass[a4paper,10pt]{article}
\usepackage[margin=1.5cm,top=2cm]{geometry}
\usepackage{graphicx}
\usepackage[hidelinks]{hyperref}
\usepackage{float}
\usepackage{parskip}

\title{\textbf{Summary -- $(ltx "$SCHEME")}}
\date{}
\begin{document}
\maketitle
HEADER

# ---------- body: one section per test ---------------------------------------

get_group() {
    case "$1" in
        toro*)       echo "toro"   ;;
        test_shear*) echo "shear"  ;;
        gresho*)     echo "gresho" ;;
        sedov*)      echo "sedov"  ;;
        *_euler)     echo "euler"  ;;
        *_ns)        echo "ns"     ;;
        *)           echo "other"  ;;
    esac
}

found=0
prev_group=""
for test in "${TESTS[@]}"; do
    outdir="$TEST_DIR/$test/outputs/$SCHEME"
    [ -d "$outdir" ] || continue

    cur_group=$(get_group "$test")

    # --- test_shear: one PDF per mesh subfolder (mesh_shear_1/2/4) ---
    if [ "$test" = "test_shear" ]; then
        shear_pdfs=()
        for mesh in mesh_shear_1 mesh_shear_2 mesh_shear_4; do
            pdf=$(find "$outdir/$mesh" -maxdepth 1 -name "*.pdf" -size +0c 2>/dev/null | head -1 || true)
            [ -n "$pdf" ] && shear_pdfs+=("$pdf")
        done
        [ "${#shear_pdfs[@]}" -gt 0 ] || continue

        [ -n "$prev_group" ] && [ "$cur_group" != "$prev_group" ] && printf '\\clearpage\n' >> "$TEX"
        found=1; prev_group="$cur_group"

        printf '\n\\section*{%s}\n' "$(ltx "$test")" >> "$TEX"
        printf '\\begin{figure}[H]\\centering\n' >> "$TEX"
        for pdf in "${shear_pdfs[@]}"; do
            local="$test-$(basename "$pdf")"
            cp "$pdf" "$WORK/$local"
            caption=$(basename "$(dirname "$pdf")" | tr '_' ' ')
            minipage "0.48" "$local" "$caption"
        done
        printf '\\end{figure}\n' >> "$TEX"
        continue
    fi

    # --- NS cylinders: images then Cp/St explicitly ---
    if [[ "$test" == *_ns ]]; then
        mapfile -t imgs < <(find "$outdir/images" -maxdepth 1 -name "*.png" -size +0c 2>/dev/null | sort)
        cp_pdf="$outdir/plot_cp.pdf"
        st_pdf="$outdir/plot_st.pdf"
        [ "${#imgs[@]}" -gt 0 ] || [ -f "$cp_pdf" ] || [ -f "$st_pdf" ] || continue

        [ -n "$prev_group" ] && [ "$cur_group" != "$prev_group" ] && printf '\\clearpage\n' >> "$TEX"
        found=1; prev_group="$cur_group"

        printf '\n\\section*{%s}\n' "$(ltx "$test")" >> "$TEX"

        if [ "${#imgs[@]}" -gt 0 ]; then
            printf '\\begin{figure}[H]\\centering\n' >> "$TEX"
            for img in "${imgs[@]}"; do
                local="$test-$(basename "$img")"
                cp "$img" "$WORK/$local"
                caption=$(basename "$img" .png | sed 's/output_-1_//; s/_/ /g')
                minipage "0.48" "$local" "$caption"
            done
            printf '\\end{figure}\n' >> "$TEX"
        fi

        if [ -f "$cp_pdf" ] || [ -f "$st_pdf" ]; then
            printf '\\begin{figure}[H]\\centering\n' >> "$TEX"
            if [ -f "$cp_pdf" ]; then
                local="$test-plot_cp.pdf"
                cp "$cp_pdf" "$WORK/$local"
                minipage "0.48" "$local" "pressure coefficient"
            fi
            if [ -f "$st_pdf" ]; then
                local="$test-plot_st.pdf"
                cp "$st_pdf" "$WORK/$local"
                minipage "0.48" "$local" "heat flux (Stanton)"
            fi
            printf '\\end{figure}\n' >> "$TEX"
        fi
        continue
    fi

    # --- toro: single multiplot PDF, shown large ---
    if [ "$cur_group" = "toro" ]; then
        mapfile -t pdfs < <(find "$outdir" -maxdepth 1 -name "*.pdf" -size +0c | sort)
        [ "${#pdfs[@]}" -gt 0 ] || continue

        [ -n "$prev_group" ] && [ "$cur_group" != "$prev_group" ] && printf '\\clearpage\n' >> "$TEX"
        found=1; prev_group="$cur_group"

        printf '\n\\section*{%s}\n' "$(ltx "$test")" >> "$TEX"
        printf '\\begin{figure}[H]\\centering\n' >> "$TEX"
        for pdf in "${pdfs[@]}"; do
            local="$test-$(basename "$pdf")"
            cp "$pdf" "$WORK/$local"
            minipage "0.70" "$local" ""
        done
        printf '\\end{figure}\n' >> "$TEX"
        continue
    fi

    # --- gresho: velocity image + radial velocity profile side by side ---
    if [ "$cur_group" = "gresho" ]; then
        velocity_img=$(find "$outdir/images" -maxdepth 1 -name "*velocity*.png" -size +0c 2>/dev/null | head -1 || true)
        if [ -f "$outdir/velocity_profile.png" ]; then profile_file="$outdir/velocity_profile.png"
        else profile_file=""; fi
        [ -n "$velocity_img" ] || [ -n "$profile_file" ] || continue

        [ -n "$prev_group" ] && [ "$cur_group" != "$prev_group" ] && printf '\\clearpage\n' >> "$TEX"
        found=1; prev_group="$cur_group"

        printf '\n\\section*{%s}\n' "$(ltx "$test")" >> "$TEX"
        printf '\\begin{figure}[H]\\centering\n' >> "$TEX"
        if [ -n "$velocity_img" ]; then
            local="$test-$(basename "$velocity_img")"
            cp "$velocity_img" "$WORK/$local"
            minipage "0.48" "$local" "velocity field"
        fi
        if [ -n "$profile_file" ]; then
            local="$test-velocity_profile.${profile_file##*.}"
            cp "$profile_file" "$WORK/$local"
            minipage "0.48" "$local" "radial velocity profile"
        fi
        printf '\\end{figure}\n' >> "$TEX"
        continue
    fi

    # --- sedov: density image + radial profile side by side ---
    if [ "$cur_group" = "sedov" ]; then
        density_img=$(find "$outdir/images" -maxdepth 1 -name "*density*.png" -size +0c 2>/dev/null | head -1 || true)
        if [ -f "$outdir/density_profile.png" ]; then profile_file="$outdir/density_profile.png"
        else profile_file=""; fi
        [ -n "$density_img" ] || [ -n "$profile_file" ] || continue

        [ -n "$prev_group" ] && [ "$cur_group" != "$prev_group" ] && printf '\\clearpage\n' >> "$TEX"
        found=1; prev_group="$cur_group"

        printf '\n\\section*{%s}\n' "$(ltx "$test")" >> "$TEX"
        printf '\\begin{figure}[H]\\centering\n' >> "$TEX"
        if [ -n "$density_img" ]; then
            local="$test-$(basename "$density_img")"
            cp "$density_img" "$WORK/$local"
            minipage "0.48" "$local" "density field"
        fi
        if [ -n "$profile_file" ]; then
            local="$test-density_profile.${profile_file##*.}"
            cp "$profile_file" "$WORK/$local"
            minipage "0.48" "$local" "radial density profile"
        fi
        printf '\\end{figure}\n' >> "$TEX"
        continue
    fi

    # --- all other tests ---
    mapfile -t pdfs < <(find "$outdir" -maxdepth 1 -name "*.pdf" -size +0c | sort)
    mapfile -t imgs < <(find "$outdir/images" -maxdepth 1 -name "*.png" -size +0c 2>/dev/null | sort)
    [ "${#pdfs[@]}" -gt 0 ] || [ "${#imgs[@]}" -gt 0 ] || continue

    [ -n "$prev_group" ] && [ "$cur_group" != "$prev_group" ] && printf '\\clearpage\n' >> "$TEX"
    found=1; prev_group="$cur_group"

    printf '\n\\section*{%s}\n' "$(ltx "$test")" >> "$TEX"

    # ParaView images first (half-page each, 2 per row)
    if [ "${#imgs[@]}" -gt 0 ]; then
        printf '\\begin{figure}[H]\\centering\n' >> "$TEX"
        for img in "${imgs[@]}"; do
            local="$test-$(basename "$img")"
            cp "$img" "$WORK/$local"
            caption=$(basename "$img" .png | sed 's/output_-1_//; s/_/ /g')
            minipage "0.48" "$local" "$caption"
        done
        printf '\\end{figure}\n' >> "$TEX"
    fi

    # gnuplot PDFs (width depends on count, larger for standalone)
    if [ "${#pdfs[@]}" -gt 0 ]; then
        n="${#pdfs[@]}"
        if   [ "$n" -le 1 ]; then w="0.45"
        elif [ "$n" -le 2 ]; then w="0.47"
        else                       w="0.30"; fi

        printf '\\begin{figure}[H]\\centering\n' >> "$TEX"
        for pdf in "${pdfs[@]}"; do
            local="$test-$(basename "$pdf")"
            cp "$pdf" "$WORK/$local"
            caption=$(basename "$pdf" .pdf | tr '_' ' ')
            minipage "$w" "$local" "$caption"
        done
        printf '\\end{figure}\n' >> "$TEX"
    fi
done

printf '\\end{document}\n' >> "$TEX"

# ---------- compile ----------------------------------------------------------

if [ "$found" -eq 0 ]; then
    echo "No outputs found for scheme '$SCHEME' in: ${TESTS[*]}" >&2
    exit 1
fi

pdflatex -interaction=nonstopmode -output-directory="$WORK" "$TEX" \
    > "$WORK/pdflatex.log" 2>&1 || {
    echo "pdflatex failed:" >&2
    tail -30 "$WORK/pdflatex.log" >&2
    exit 1
}

cp "$WORK/summary.pdf" "$OUT_PDF"
echo "Generated: $OUT_PDF"
