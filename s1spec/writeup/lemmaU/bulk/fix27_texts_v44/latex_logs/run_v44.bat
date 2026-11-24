@echo off
set PATH=C:\Users\micahw\AppData\Local\Programs\MiKTeX\miktex\bin\x64;%PATH%
set F=%1
pdflatex -interaction=nonstopmode -halt-on-error %F%.tex > nul
findstr /c:"\bibdata" %F%.aux > nul && bibtex %F% > nul
pdflatex -interaction=nonstopmode -halt-on-error %F%.tex > nul
pdflatex -interaction=nonstopmode -halt-on-error %F%.tex > nul
echo rc=%ERRORLEVEL%
