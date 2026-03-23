#!/bin/bash
module load ufrc mkl/2023.2.0 gcc/12.2.0 openmpi/4.1.5 python/3.11 cmake/3.26.4
k=10
for m in {10..2000..10}
do
    let n=$m/$k
    var1=$(printf "\'PFData%03d\'" "$n")
    mpirun /blue/michael.tonks/share/trgapp/gatorapp-opt -i PolycrystalGrainGrowth.i file=$var1 randseed=$m

    
done