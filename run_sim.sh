#!/bin/bash

make clean
make dust_simulation
((X = 5000))
((Y = 5000))
((num = 2500))
((min = 1))
((max = 30))
((t = 3600))
((stick = 1))
((splt = 0))
((mrg = 0))
((spacing = 6))
((gap = 9))
((depth = 3))
((vx = 4*$max))
((vy = 4*$max))
echo "Launching dust filtration simulation(s)..."
for ((id=0; id<1; id=id+1))
do
	echo "Running ./dust_simulation "$X" "$Y" "$vx" "$vy" "$num" "$min" "$max" "$t" "$stick" "$splt" "$mrg" "$spacing" "$gap" "$depth" "$id
	./dust_simulation $X $Y $vx $vy $num $min $max $t $stick $splt $mrg $spacing $gap $depth $id
done
echo "All Done!"
