Implementation of [rheotool](https://github.com/fppimenta/rheoTool).

Configure ParaFoam with arbitrary version of Paraview:
* Modify $etc/config.sh/paraview and $etc/config.csh/paraview with the correct Paraview version and path.
* ParaFoam command is tricky, but not necessary. Paraview is faster indicated [here](http://www.wolfdynamics.com/training/introOF8/supplement_paraview_parafoam.pdf).

Rheotool installation in OpenFOAM 9:
* The environmental variables for PETSC should be added in $etc/bashrc
