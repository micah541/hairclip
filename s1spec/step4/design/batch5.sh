cd /home/claude/s1spec/step4/design
export N0=32768
LEVELS=100,200,400,800,1600 python3 run3.py clip_24_0.4_0.05_0.85 clip:24:0.4:0.05:0.85 1700
CDT=5e-3 TOL=1e-12 LEVELS=400,800 python3 run3.py fine_clip_48_0.4_0.05 clip:48:0.4:0.05 900
