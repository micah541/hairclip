cd /home/claude/s1spec/step4/design
export LEVELS=100,200,400,800
for spec in clip:8:0.3:0.05 clip:5:0.3:0.05 clip:3:0.3:0.05 clip:12:0.2:0.05 clip:12:0.4:0.05 clip:12:0.3:0.1 clip:12:0.3:0.05:1.5; do
  tag=$(echo $spec | tr ':' '_'); python3 run3.py $tag $spec 900
done
