cd /home/claude/s1spec/step4/design
export LEVELS=100,200,400,800
python3 run3.py disk_R2.5 disk:2.5:0.05:0.3:0 900
python3 run3.py clip_a12_c0.3 clip:12:0.3:0.05 900
