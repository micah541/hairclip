@echo off
set F=%1
pdflatex -interaction=nonstopmode -halt-on-error "\input{%F%}" > nul
findstr /c:"\bibdata" %F%.aux > nul && bibtex %F% > nul
pdflatex -interaction=nonstopmode -halt-on-error "\input{%F%}" > nul
pdflatex -interaction=nonstopmode -halt-on-error "\input{%F%}" > nul
