#!/bin/bash

nc=100
ni=20

# loop through each center/image
for i in `seq $ni`
do

x=`head -$i centers.txt | tail -1 | awk '{print $2}'`
y=`head -$i centers.txt | tail -1 | awk '{print $3}'`

# paramter file
cat > toy.prm << EOF
begin

timestep 1e-6
runtime	2
damping 1
strength 0.0313716
initial $x $y
report 0.001

US 0
center $x $y

end

EOF

# run simulations: loop through each copy (for each image)
for j in `seq $nc`
do
    ./toy_md "out/$i.$j.gz" $((i*nc+j)) < toy.prm &
done

while pgrep toy_md > /dev/null; do sleep 1; done

done

# gunzip -c out/1.2.gz | head
