#!/bin/bash
# for i in $(foamListTimes -case processor0); do
#     reconstructPar -time ${i}
#     rm -r processo*/${i}
# done
#
# rm -rf processor*

for i in */; do zip -0 -r "${i%/}.zip" "$i" & done; wait
