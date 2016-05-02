#!/bin/sh
make clean
make dust_simulation

echo "Doing dust simulation..."

./dust_simulation 2000 2000 80 80 2000 10 40 3500 1 0 0 10 10 10 0

echo "All Done!"
