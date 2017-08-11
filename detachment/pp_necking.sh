#!/bin/bash

## This script merges the particle output data for each timestep
## and then computes the minimum distance between the two particle
## at the same depth.

## User-defined parameters
## The directory into which ASPECT output is written
dir=detachment

## Remove previously created postprocessing files
rm $dir/particles/necking.txt

## Create output file and write header
#echo "#time [s], Output_step [-], necking width [m], necking depth [m]" > $dir/particles/necking.txt

## Loop over all output times
for t in {00..25}

do

 ## Remove any previously created postprocessing files
 rm $dir/particles/particle-000$t.txt

 ## Append the entries for each processor to one file 
 grep -h -v "#" $dir/particles/particle-000$t.00*.txt >> $dir/particles/particle-000$t.txt

 ## Sort according to depth of the particles 
 sort -k 2n $dir/particles/particle-000$t.txt > $dir/particles/particle-000$t.sorted.txt

 ## Write the particles with the same depth on one line
 paste -s -d '\t\n' $dir/particles/particle-000$t.sorted.txt > $dir/particles/particle-000$t.sorted.oneline.txt

 ## Compute which two particles are closest, and output their distance and depth
 awk 'function abs(v){return v<0. ? -v : v}BEGIN{dist=1e10;D=660000.;depth=1e10}{if(abs($4-$1)<dist) {dist=abs($4-$1); depth=D-$2}}END{print '$t',dist,depth}' $dir/particles/particle-000$t.sorted.oneline.txt >> $dir/particles/necking.txt

done

## Add in the model time 
## NB: ASPECT 1.5 time in statistics file is off by one timestep,
## so we use the output times of the visualization file
grep timestep $dir/solution.pvd | awk 'BEGIN{FS="\""}{print $2}' | paste -d ' ' - $dir/particles/necking.txt > $dir/particles/necking_evolution.txt
