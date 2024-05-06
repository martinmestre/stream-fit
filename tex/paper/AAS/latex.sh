#!/bin/bash
#
#


pdflatex sample631.tex;
bibtex sample631.aux;
pdflatex sample631.tex;
pdflatex sample631.tex;




