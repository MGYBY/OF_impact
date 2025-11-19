# restore0Dir
cp -r 0_org 0

fluent3DMeshToFoam ansys_mesh.msh
# blockMesh
createPatch -overwrite
# exit

# renumberMesh -overwrite

# I.C.
funkySetFields -time 0

# exit

# # # # # for parallel running:
decomposePar

# exit

postProcess

# mpirun --oversubscribe -np 15 rheoInterFoam -parallel
mpirun --oversubscribe -np 15 interIsoDeboRheoFoam -parallel
# rheoInterFoam
# interIsoDeboRheoFoam
#------------------------------------------------------------------------------
