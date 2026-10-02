#!/bin/bash
# Combines gresho_quad/gresho_tri/convergence_gresho_{quad,tri} WIP images into
# one summary.pdf. Needs pdflatex, which is NOT available on Curta — run this
# locally, after rsync'ing the cluster outputs/ back (see ../README.md).
set -e

ROOT=$(cd "$(dirname "$0")"; pwd)
CASEROOT=$(dirname "$ROOT")
OUTDIR="$ROOT/outputs/WIP"
mkdir -p "$OUTDIR"

Q="$CASEROOT/gresho_quad/outputs/WIP/images"
T="$CASEROOT/gresho_tri/outputs/WIP/images"
CQ="$CASEROOT/convergence_gresho_quad/outputs/WIP/images"
CT="$CASEROOT/convergence_gresho_tri/outputs/WIP/images"

TEX="$OUTDIR/summary.tex"

cat > "$TEX" << TEX_EOF
\documentclass[a4paper]{article}
\usepackage[margin=0.5cm]{geometry}
\usepackage{graphicx}
\usepackage[T1]{fontenc}

\newcommand{\safeimg}[2]{%
  \IfFileExists{#1}%
    {\includegraphics[width=#2,height=#2,keepaspectratio]{#1}}%
    {\fbox{\parbox{#2}{\centering\small N/A}}}%
}

\begin{document}
\pagestyle{empty}

\begin{center}{\large\textbf{WIP}}\end{center}
\vspace{4pt}

\begin{center}\textbf{Quad mesh}\end{center}\vspace{2pt}
\noindent
\begin{tabular}{@{}ccc@{}}
  \safeimg{$Q/gresho_quad_density.png}{0.31\textwidth} &
  \safeimg{$Q/gresho_quad_pressure.png}{0.31\textwidth} &
  \safeimg{$Q/gresho_quad_velocity.png}{0.31\textwidth} \\\\[4pt]
  \safeimg{$Q/velocity_profile.pdf}{0.31\textwidth} &
  \safeimg{$Q/residual.pdf}{0.31\textwidth} &
  \safeimg{$CQ/convergence.pdf}{0.31\textwidth} \\\\
\end{tabular}

\vspace{8pt}
\begin{center}\textbf{Tri mesh}\end{center}\vspace{2pt}
\noindent
\begin{tabular}{@{}ccc@{}}
  \safeimg{$T/gresho_tri_density.png}{0.31\textwidth} &
  \safeimg{$T/gresho_tri_pressure.png}{0.31\textwidth} &
  \safeimg{$T/gresho_tri_velocity.png}{0.31\textwidth} \\\\[4pt]
  \safeimg{$T/velocity_profile.pdf}{0.31\textwidth} &
  \safeimg{$T/residual.pdf}{0.31\textwidth} &
  \safeimg{$CT/convergence.pdf}{0.31\textwidth} \\\\
\end{tabular}

\end{document}
TEX_EOF

cd "$OUTDIR"
pdflatex -interaction=batchmode summary.tex > pdflatex.log 2>&1

date > summary.timestamp
echo "-> $OUTDIR/summary.pdf"
