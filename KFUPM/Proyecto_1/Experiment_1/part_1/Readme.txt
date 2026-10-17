========================================================================
VERTEX-BASED STATIONARY LINEAR ELASTICITY
========================================================================

Developed at King Fahd University of Petroleum and Minerals (KFUPM)
Dhahran, Saudi Arabia -- IRC for Industrial Nuclear Energy

------------------------------------------------------------------------
DESCRIPTION
------------------------------------------------------------------------
This MATLAB project compares a vertex-based discretization of stationary
linear elasticity using random and regular interior meshes of the unit disk.
Both families share equally spaced circular boundary points. The regular
interior uses a Cartesian lattice with consistent triangulation diagonals;
its actual node count may differ slightly from the requested count. The
computational domain is the inscribed boundary polygon.

------------------------------------------------------------------------
EXECUTION
------------------------------------------------------------------------
Compatible with MATLAB R2017b. Open the project directory in MATLAB and run:
    main

When prompted, enter a target node count $N \geq 30$. A scalar $N$ generates 
up to three refinement levels, with repeated counts removed.

------------------------------------------------------------------------
MATHEMATICAL SETUP
------------------------------------------------------------------------
The exact displacement is $u(x,y) = [x^2+y^2; 2x y]^T$. With Lamé parameters 
$\lambda = 2$ and $\mu = 1$, the body force is $f = [-16; 0]^T$. Dirichlet data 
use the exact displacement; Neumann traction is $\sigma_{\text{exact}} n_{\text{edge}}$ 
on the upper polygon edges.

------------------------------------------------------------------------
RESULTS AND DIAGNOSTICS
------------------------------------------------------------------------
The program compares analytical and discrete fields and reports displacement,
tensor, and gradient errors, including diagonal and off-diagonal components.
Tensor errors use the full Frobenius norm; weighted $L_2$ errors use nodal
quadrature. Observed convergence orders are computed separately for each mesh
family using the maximum actual triangulation edge length. Invalid orders
are $\text{NaN}$; rates may fluctuate or be negative.

Figures show the finest level, including a 2-by-3 comparison of discrete
deformations, undeformed disks, and analytical deformations. Colors encode
the original angle and radius; deformations use scale factor 1.
Data from every level, error tables, convergence orders, and equilibrium
residuals are saved to `Data_files.mat` in the current MATLAB folder.

Note: The free system uses `pinv(full(KFF))`. This dense pseudoinverse can 
become computationally expensive for large node counts.