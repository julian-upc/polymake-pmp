#!/bin/bash
# This script takes as input 
# . the filename of a reference polytope P that includes a GROUP property (e.g., QUOTIENT_SPACE->SYMMETRY_GROUP)
# . the dimension of P
# . a directory name inside which the results of computation will be stored.
# it iteratively constructs the equivalence classes of k-simplices constructed from the vertices of P
# such that the computations for a fixed k are executed in parallel.

REF_POLY=$1 # the reference polytope
DIM=$2      # the dimension
DATA_DIR=$3 # where to put the results

N_CORES=`cat /proc/cpuinfo | grep processor | wc -l`   # the number of available computation kernels


polymake 'my $r=load("'$REF_POLY'"); make_facet_makefile($r, $r, "'$REF_POLY'", out_dir=>"'$DATA_DIR'");'
cd $DATA_DIR/dim1_; make -j$N_CORES; cd ../..
for ((i=2; i<=$DIM; i++)) do
    polymake 'my $r=load("'$REF_POLY'"); my $c=load("'$DATA_DIR'/step'$((i-1))'.poly"); make_facet_makefile($c, $r, "'$REF_POLY'", out_dir=>"'$DATA_DIR'");'
    cd $DATA_DIR/dim${i}_; make -j$N_CORES; cd ../..
done
