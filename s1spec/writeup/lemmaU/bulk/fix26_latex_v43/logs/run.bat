@echo off
set F=%1
pdflatex -interaction=nonstopmode -halt-on-error "\providecommand{\key}[1]{\texttt{\detokenize{#1}}}\input{%F%}" > nul
findstr /c:"\bibdata" %F%.aux > nul && bibtex %F% > nul
pdflatex -interaction=nonstopmode -halt-on-error "\providecommand{\key}[1]{\texttt{\detokenize{#1}}}\input{%F%}" > nul
pdflatex -interaction=nonstopmode -halt-on-error "\providecommand{\key}[1]{\texttt{\detokenize{#1}}}\input{%F%}" > nul
