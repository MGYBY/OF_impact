//- Create segregated solvers for p_rgh|U if needed
// autoPtr<sparseSolver<scalar> > spSolverprgh
autoPtr<sparseSolver<scalar> > spSolverp
(
  // sparseSolver<scalar>::New(p_rgh, mesh, mesh.solutionDict()).ptr()
  sparseSolver<scalar>::New(p, mesh, mesh.solutionDict()).ptr()
);

autoPtr<sparseSolver<vector> > spSolverU
(
  sparseSolver<vector>::New(U, mesh, mesh.solutionDict()).ptr()
);


