#!/bin/sh
# one command: rebuild (add the word zip to also write ../kit_step4b.zip) kit_step4b from the proof tree and zip it (files at the zip root)
set -e
cd "$(dirname "$0")"
python3 build.py
cd ../kit_step4b
find . -name __pycache__ -prune -exec rm -rf {} + ; rm -rf kit_out
[ "$1" = zip ] || { echo 'built (no zip; run: sh refresh.sh zip)'; exit 0; }
find . -name __pycache__ -prune -exec rm -rf {} + ; rm -rf kit_out
rm -f ../kit_step4b.zip && python3 -c "
import zipfile,os
z=zipfile.ZipFile('../kit_step4b.zip','w',zipfile.ZIP_DEFLATED)
for r,d,f in os.walk('.'):
    for n in f: p=os.path.join(r,n); z.write(p,os.path.relpath(p,'.'))
print('zip', os.path.getsize('../kit_step4b.zip'))"
