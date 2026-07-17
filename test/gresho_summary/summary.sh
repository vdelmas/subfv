#!/bin/bash
set -e

SCHEME=$1
ROOT=$(cd "$(dirname "$0")"; pwd)
TESTROOT=$(dirname "$ROOT")
OUTDIR="$ROOT/outputs/$SCHEME"
mkdir -p "$OUTDIR"

Q="$TESTROOT/gresho_quad/outputs/$SCHEME/images"
T="$TESTROOT/gresho_tri/outputs/$SCHEME/images"
CQ="$TESTROOT/convergence_gresho_quad/outputs/$SCHEME/images"
CT="$TESTROOT/convergence_gresho_tri/outputs/$SCHEME/images"

# Escape underscores for LaTeX text mode
SCHEME_LATEX=$(echo "$SCHEME" | sed 's/_/\\_/g')

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

\begin{center}{\large\textbf{$SCHEME_LATEX}}\end{center}
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
