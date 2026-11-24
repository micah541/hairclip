for k in 4000 5000 6000 7000 8000 10000; do python3 drive.py snap_main2_$k.npz 8 main_$k > drive_main_$k.out 2>&1; done
