for pair in main_1600:snap_main_1600 main_2000:snap_main_2000 main_3000:snap_main_3000 main_4000:snap_main2_4000 main_5000:snap_main2_5000 main_6000:snap_main2_6000 main_7000:snap_main2_7000 main_8000:snap_main2_8000 main_10000:snap_main2_10000 big_8000:snap_big_8000 res20_8000:snap_res20_8000 big_10000:snap_big_10000 fine_2000:snap_fine_2000 fine_3000:snap_fine_3000 fine_4000:snap_fine_4000 hr_4000:snap_hr_4000; do
  t=${pair%%:*}; s=${pair##*:}; python3 drive.py $s.npz 16 $t > drive_$t.out 2>&1
done
