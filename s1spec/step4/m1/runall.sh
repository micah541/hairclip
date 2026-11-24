for k in 2000 3000 4000 5000 6000 7000 8000 10000; do
  f=snap_main_$k.npz; [ -f snap_main2_$k.npz ] && f=snap_main2_$k.npz
  python3 drive.py $f 16 main_$k > drive_main_$k.out 2>&1
done
python3 drive.py snap_main_1600.npz 16 main_1600 > drive_main_1600.out 2>&1
