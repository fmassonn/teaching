#!/bin/bash

set -o nounset
set -o errexit

for tag in `cat listTime.txt`
do
    echo $tag
    convert ./figs/z500_${tag}.png ./figs/t500_${tag}.png +append 1.png
    convert ./figs/z1000_${tag}.png ./figs/t1000_${tag}.png +append 2.png
    convert 1.png 2.png -append ./figs/all_${tag}.png
done

ffmpeg -framerate 30 -pattern_type glob -i './figs/all_*.png'  -c:v libx264 -pix_fmt yuv420p out.mp4
