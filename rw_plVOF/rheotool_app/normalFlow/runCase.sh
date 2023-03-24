cp -r 0_org 0
blockMesh
# setFields
funkySetFields -time 0
decomposePar
# mpirun -np 4 interCyclicFoam -parallel
# interCyclicFoam
# rheoInterFoam
mpirun -np 4 rheoInterFoam -parallel
