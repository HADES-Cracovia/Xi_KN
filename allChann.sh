#!/bin/bash

evt=4000000
arr=( 090 010 021 022 023 062 100 )
for i in ${arr[@]}; do
    echo $i
    ./vertex_rec out/pluto_chan_${i}_events_10000_seed_000_1_dst_p4500p.root -d ./ -o output_${i}.root -e $evt > /dev/null &
done
