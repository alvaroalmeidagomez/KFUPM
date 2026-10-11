% =========================================================================
% Developed at King Fahd University of Petroleum and Minerals (KFUPM)
% Dhahran, Saudi Arabia
% IRC for Industrial Nuclear Energy
% =========================================================================
% =========================================================================
% Developed at King Fahd University of Petroleum and Minerals (KFUPM)
% Dhahran, Saudi Arabia
% IRC for Industrial Nuclear Energy
% =========================================================================

% RUN_DISK_ELASTICITY_EXPERIMENT Numerical experiment for the vertex-based 
% discretization of stationary linear elasticity on the unit disk.
%
% Description:
%   This script constructs two point sets with N points, computes their 
%   Delaunay triangulations, evaluates the nodal data, assembles the 
%   vertex-based discrete weak formulation, and solves the system for both 
%   triangulations.
%
% Execution Notes:
%   - Uses MATLAB R2017b-compatible plotting syntax.
%   - Prompts the user interactively for the number of points (N >= 30).
%   - Automatically generates and displays all figures.
%   - Utilizes the dense Moore-Penrose pseudoinverse to solve the free 
%     displacement system.
%
% Mathematical Setup:
%   The exact displacement field is given by:
%       u_exact(x,y) = [x^2 + y^2; 2*x*y]
%
%   Using Lame parameters lambda = 2 and mu = 1, the corresponding source 
%   term (f) and traction boundary condition (t) are:
%       f = [-16; 0]
%       t = [12*x^2 + 4*y^2; 16*x*y]

% Deformation figures display Phi(X) = X + U(X), with scale equal to 1.
% One large 2-by-3 figure groups all six deformation-related images.
% Rows: random and quasi-uniform. Columns: discrete deformation,
% original disk, and analytical deformation. Reference-node colors
% are preserved throughout. These six panels show only points.
% Hue encodes the original polar angle; saturation and brightness encode
% the original radius. A two-dimensional color key explains both variables.
% Both point sets use the SAME equally spaced boundary; only their
% interior distributions differ. Boundary nodes are drawn last.
%
% Gradient diagnostics use G(i,j,a) = partial u_i / partial x_j at node a.
% Nodal full/diagonal/off-diagonal errors are Euclidean norms of the
% corresponding entries. Their means are arithmetic means over all nodes.

    clc;
    close all;
    clear all
    
    run_disk_elasticity_experiment