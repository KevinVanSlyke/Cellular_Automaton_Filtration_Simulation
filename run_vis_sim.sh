#!/bin/bash

make clean
make visual_dust_simulation
((X = 1000))
((Y = 1000))
((num = 100))
((min = 1))
((max = 10))
((t = 10000))
((stick = 1))
((splt = 1))
((mrg = 1))
((spacing = 6))
((gap = 9))
((depth = 3))
((Px = 4*$max))
((Py = 4*$max))
((nPy = 4*$max))
echo "Launching dust filtration simulation(s)..."
for ((id=0; id<1; id=id+1))
do
	echo "Running ./visual_dust_simulation "$X" "$Y" "$Px" "$Py" "$nPy" "$num" "$min" "$max" "$t" "$stick" "$splt" "$mrg" "$spacing" "$gap" "$depth" "$id
	./visual_dust_simulation $X $Y $Px $Py %nPy $num $min $max $t $stick $splt $mrg $spacing $gap $depth $id
done
echo "All Done!"
