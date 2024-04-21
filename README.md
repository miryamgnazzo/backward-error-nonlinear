# backward-error-nonlinear
This <code>backward-error-nonlinear</code> repository contains the MATLAB code for the computation of the backward error
associated with a set of eigenpairs of a nonlinear eigenvalue problem. The computation of both the structured and
unstructured backward error follows the methodology in [1].

## Main functions:
* **be_unstructured.m**: computes the unstructured backward error for a set of eigenpairs;
* **be_unstructured_bound.m**: computes the bounds for the unstructured backward error;
* **be_linear_structured.m**: computes the structured backward error for a prescribed linear structure;
* **be_linear_structured_bound.m**: computes the bounds for the linear structured backward error;
* **be_symmetric.m**: computes the symmetric backward error;
* **be_symmetric_bound.m**: computes the bounds for the symmetric backward error; 
* **be_riemannian.m**: computes an upper bound for the structured backward error,
                       for nonlinear structures (sparsity structures, low-rank, multiple of the identity),
                        by a Riemannian optimization based method (Section 3.3, [1]).

The folder <code>examples/</code> collects the codes for the numerical experiments provided in Section 4 [1].

[1] M. Gnazzo, L. Robol, Backward errors for multiple eigenpairs in structured and unstructured nonlinear eigenvalue problems,
available on arXiv (2024).

## Dependencies
* The use of the Riemannian optimization-based approach requires the [manopt package](https://www.manopt.org/index.html).
* A few nonlinear eignvalue problems are taken from the collection [nlevp](https://github.com/ftisseur/nlevp?tab=readme-ov-file).
