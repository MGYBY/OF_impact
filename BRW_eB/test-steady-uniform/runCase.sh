cp -r 0_org 0
blockMesh
# setFields
# funkySetFields -time 0
setFields
decomposePar
# mpirun -np 4 interCyclicFoam -parallel
# interCyclicFoam
# rheoInterFoam
# mpirun -np 6 interFoam -parallel
mpirun -np 6 foamRun -solver incompressibleVoF -parallel
