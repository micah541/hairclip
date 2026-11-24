@echo off
set F=%1
C:\Users\micahw\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe -interaction=nonstopmode -halt-on-error "\input{%F%}" > nul
findstr /c:"\bibdata" %F%.aux > nul && C:\Users\micahw\AppData\Local\Programs\MiKTeX\miktex\bin\x64\bibtex.exe %F% > nul
C:\Users\micahw\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe -interaction=nonstopmode -halt-on-error "\input{%F%}" > nul
C:\Users\micahw\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe -interaction=nonstopmode -halt-on-error "\input{%F%}" > nul
