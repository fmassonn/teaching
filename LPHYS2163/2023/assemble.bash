#!/bin/bash

set -o nounset
set -o errexit

for tag in `cat listTime.txt`
do
    convert ./figs/z500_${tag}.png ./figs/t500_${tag}.png +append ./figs/all_${tag}.png
done
