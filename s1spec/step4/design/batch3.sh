cd /home/claude/s1spec/step4/design
export LEVELS=100,200,400,800
for spec in clip:12:0.5:0.05 clip:12:0.6:0.05 clip:12:0.3:0.05:0.7 clip:12:0.4:0.03; do
  tag=$(echo $spec | tr ':' '_'); python3 run3.py $tag $spec 900
done
export N0=32768
for spec in clip:36:0.4:0.05 clip:48:0.4:0.05; do
  tag=$(echo $spec | tr ':' '_'); python3 run3.py $tag $spec 900
done
