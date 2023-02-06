cp -r 0_org 0
blockMesh
setFields
decomposePar
mpirun -np 4 interCyclicFoam -parallel
