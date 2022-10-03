3D impact problem.
1. Turbulent model: k-o SST
2. Half-domain
3. Circular cylinder. Special attention to meshing.
4. interIsoFoam
5. RDF_PLIC
6. swak4Foam for IC and BC from 1D simulations.
7. enable hyprethread for parallelism:
```
mpirun --oversubscribe -np 12 interFoam -parallel
```
