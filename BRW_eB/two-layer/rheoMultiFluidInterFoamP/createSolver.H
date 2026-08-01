// Create segregated sparse solvers for physical pressure p and U.
autoPtr<sparseSolver<scalar> > spSolverp
(
    sparseSolver<scalar>::New(p, mesh, mesh.solutionDict()).ptr()
);

autoPtr<sparseSolver<vector> > spSolverU
(
    sparseSolver<vector>::New(U, mesh, mesh.solutionDict()).ptr()
);
