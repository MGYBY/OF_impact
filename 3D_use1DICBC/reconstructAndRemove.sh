#!/bin/bash
for i in $(foamListTimes -case processor0); do
    reconstructPar -time ${i}
    rm -r processo*/${i}
done
