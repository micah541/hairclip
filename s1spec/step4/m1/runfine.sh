for k in 2000 3000 4000; do python3 drive.py snap_fine_$k.npz 16 fine_$k > drive_fine_$k.out 2>&1; done
