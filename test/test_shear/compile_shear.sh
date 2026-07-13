#!/bin/bash
set -e
pdflatex -interaction=nonstopmode -halt-on-error shear_comparison.tex
pdflatex -interaction=nonstopmode -halt-on-error shear_comparison.tex
