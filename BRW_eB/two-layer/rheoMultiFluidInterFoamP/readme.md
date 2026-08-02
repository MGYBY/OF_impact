`rheoMultiFluidInterFoam` modified to cyclic BCs.

Sometimes needs the following variables:
```
source /opt/openfoam7/etc/bashrc

export FOAM_USER_SRC="$WM_PROJECT_USER_DIR/src"
export FOAM_USER_LIBRARIES="$FOAM_USER_SRC/libs"
export FOAM_USER_SOLVERS="$FOAM_USER_SRC/solvers"
export FOAM_USER_THIRDPARTY="$WM_PROJECT_USER_DIR/ThirdParty"

export EIGEN_RHEO="$FOAM_USER_THIRDPARTY/Eigen3.2.9"
export PETSC_DIR="$FOAM_USER_THIRDPARTY/petsc-3.15.0"
export PETSC_ARCH=arch-linux-c-opt

export LD_LIBRARY_PATH="\
$PETSC_DIR/$PETSC_ARCH/lib:\
$FOAM_USER_LIBBIN:\
${LD_LIBRARY_PATH:-}"
```

The `HerschelBulkley.C` file is in `\src\libs\constitutiveEquations\constitutiveEqs\GNF\HerschelBulkley\` rheoTool path. Replacement commands:
```
RHEO_LIBS="$(
    readlink -f \
    "$FOAM_USER_SRC/libs/libRheoTool"
)"
cd "$RHEO_LIBS/constitutiveEquations"
wclean
wmake libso
```
