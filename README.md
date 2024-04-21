# backward-error-nonlinear
This <code>backward-error-nonlinear</code> repository contains the MATLAB code for the computation of the backward error
associated with a set of eigenpairs of a nonlinear eigenvalue problem. The computation of both the structured and
unstructured backward error follows the methodology in [1].

## Main functions:
* **be_unstructure.m**



[1] M. Gnazzo, L. Robol, Backward errors for multiple eigenpairs in structured and unstructured nonlinear eigenvalue problems,
available on arXiv (2024).

The use of the Riemannian optimization-based approach requires the [manopt package](https://www.manopt.org/index.html).

The folder <code>examples/</code> collects a few examples.




