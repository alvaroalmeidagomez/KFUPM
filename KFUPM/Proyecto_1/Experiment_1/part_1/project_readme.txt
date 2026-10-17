=========================================================================
Developed at King Fahd University of Petroleum and Minerals (KFUPM)
Dhahran, Saudi Arabia
IRC for Industrial Nuclear Energy
=========================================================================

RUN_DISK_ELASTICITY_EXPERIMENT
Numerical experiment for the vertex-based discretization of stationary linear elasticity on the unit disk.

-------------------------------------------------------------------------
Description:
-------------------------------------------------------------------------
This project constructs two point sets with N points, computes their 
Delaunay triangulations, evaluates the nodal data, assembles the 
vertex-based discrete weak formulation, and solves the system for both 
triangulations.

-------------------------------------------------------------------------
Execution Instructions:
-------------------------------------------------------------------------
1. Ensure you have MATLAB R2017b (or a compatible version) installed.
2. Open MATLAB and navigate to the project directory.
3. Run the main script from the MATLAB command window:
   >> main
4. The program will prompt you interactively to enter the number of points (N >= 30).
5. Upon execution, all deformation figures and diagnostic plots are automatically generated and displayed.
6. After execution completes, a file named `data_files.m` is saved in the directory, containing all the generated experiment information and data variables.

-------------------------------------------------------------------------
Execution Notes:
-------------------------------------------------------------------------
- Uses MATLAB R2017b-compatible plotting syntax.
- Utilizes the dense Moore-Penrose pseudoinverse to solve the free displacement system.

-------------------------------------------------------------------------
Mathematical Setup:
-------------------------------------------------------------------------
The exact displacement field is given by:
    u_exact(x,y) = [x^2 + y^2; 2*x*y]

Using Lame parameters lambda = 2 and mu = 1, the corresponding source 
term (f) and traction boundary condition (t) are:
    f = [-16; 0]
    t = [12*x^2 + 4*y^2; 16*x*y]

-------------------------------------------------------------------------
Additional Details:
-------------------------------------------------------------------------
- Deformation figures display Phi(X) = X + U(X), with scale equal to 1.
- One large 2-by-3 figure groups all six deformation-related images (Rows: random and quasi-uniform; Columns: discrete deformation, original disk, and analytical deformation).
- Hue encodes the original polar angle; saturation and brightness encode the original radius.
- Gradient diagnostics use G(i,j,a) = partial u_i / partial x_j at node a. Nodal full/diagonal/off-diagonal errors are Euclidean norms of the corresponding entries.