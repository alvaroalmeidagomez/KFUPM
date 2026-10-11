function results = run_disk_elasticity_experiment
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
    
disp('=========================================================================');
disp('Developed at King Fahd University of Petroleum and Minerals (KFUPM)');
disp('Dhahran, Saudi Arabia');
disp('IRC for Industrial Nuclear Energy');
disp('=========================================================================');

%% Parameters
while true

    textN = input('Enter the number of points N (integer >= 30): ', 's');
    N = str2double(strtrim(textN));

    if ~isnan(N) && isfinite(N) && N == fix(N) && N >= 30
        break;
    end

    fprintf('Invalid input. N must be an integer >= 30 (at least 3 boundary points).\n\n');
end

    makePlots = true;
    useMoorePenrose = true;

    % Approximately 10% of the nodes lie on the boundary. The rounded
    % count is at least 3 because N >= 30. N itself is not changed.
    Nboundary = round(0.10*N);
    lambda = 2;
    mu = 1;

    rng(1,'twister');

    %% Common EQUALLY SPACED boundary points on the unit circle
    % Both distributions receive these exact same coordinates, with
    % angular spacing 2*pi/Nboundary and equal arc and chord lengths.
    % Only the interior nodes of P_rand are sampled randomly.
    thetaBoundary = (0:Nboundary-1)'*(2*pi/Nboundary);
    boundaryPoints = [cos(thetaBoundary), sin(thetaBoundary)];

    %% Point distributions
    P_rand = make_point_set(N, 'random', boundaryPoints);
    P_eq   = make_point_set(N, 'quasiuniform', boundaryPoints);

    %% Delaunay triangulations and geometric weights
    meshRand = prepare_mesh(P_rand);
    meshEq   = prepare_mesh(P_eq);

    %% Discrete weak solutions
    solRand = solve_discrete_problem(meshRand, lambda, mu);
    solEq   = solve_discrete_problem(meshEq,   lambda, mu);

    %% Store results
    results.parameters.N = N;
    results.parameters.Nboundary = Nboundary;
    results.parameters.boundaryFraction = Nboundary/N;
    results.parameters.boundaryPoints = boundaryPoints;
    results.parameters.boundaryAngleStep = 2*pi/Nboundary;
    results.parameters.lambda = lambda;
    results.parameters.mu = mu;
    results.parameters.makePlots = makePlots;
    results.parameters.useMoorePenrose = useMoorePenrose;
    results.parameters.deformationScale = 1;
    results.random.mesh = meshRand;
    results.random.solution = solRand;
    results.quasiuniform.mesh = meshEq;
    results.quasiuniform.solution = solEq;

    %% Numerical summary
    fprintf('\nVertex-based elasticity experiment on the unit disk\n');
    fprintf(['N = %d, boundary vertices = %d (%.2f%%), ' ...
             'lambda = %.1f, mu = %.1f\n\n'], ...
        N, Nboundary, 100*Nboundary/N, lambda, mu);

    print_summary('Random', meshRand, solRand);
    print_summary('Quasi-uniform', meshEq, solEq);

    %% Figures (always displayed)
    plot_meshes(meshRand, meshEq);
    plot_solutions(meshRand, solRand, meshEq, solEq);
    plot_deformations(meshRand, solRand, meshEq, solEq);
    plot_gradient_errors(meshRand, solRand, meshEq, solEq);
    
    
    
    save('Data_files.mat');
end


function P = make_point_set(N, distribution, boundaryPoints)
%MAKE_POINT_SET Create N points with a shared polygonal boundary.

    Nboundary = size(boundaryPoints,1);
    Ninterior = N - Nboundary;

    % Sampling within the polygon's incircle keeps all interior points
    % strictly inside the regular polygon inscribed in the unit circle.
    rmax = cos(pi/Nboundary)*(1 - 1e-12);

    if strcmpi(distribution,'random')
        interior = zeros(Ninterior,2);
        nfilled = 0;
        while nfilled < Ninterior
            nremaining = Ninterior - nfilled;
            batchSize = max(100,2*nremaining);
            candidate = 2*rand(batchSize,2) - 1;
            inside = sum(candidate.^2,2) <= rmax^2;
            candidate = candidate(inside,:);

            nnew = min(size(candidate,1),nremaining);
            if nnew > 0
                interior(nfilled+1:nfilled+nnew,:) = candidate(1:nnew,:);
                nfilled = nfilled + nnew;
            end
        end

    elseif strcmpi(distribution,'quasiuniform') || strcmpi(distribution,'eq')
        % Fermat spiral with approximately uniform area density.
        k = (1:Ninterior)';
        goldenAngle = pi*(3 - sqrt(5));
        theta = goldenAngle*(k-1);
        radius = rmax*sqrt((k - 0.5)/Ninterior);
        interior = [radius.*cos(theta), radius.*sin(theta)];

    else
        error('Unknown point distribution: %s',distribution);
    end

    P = [boundaryPoints; interior];
end


function mesh = prepare_mesh(P)
%PREPARE_MESH Construct the mesh, boundary partition, and vertex weights.

    n = size(P,1);
    DT = delaunayTriangulation(P);
    T = DT.ConnectivityList;
    boundaryEdges = freeBoundary(DT);

    isBoundary = false(n,1);
    isBoundary(unique(boundaryEdges(:))) = true;

    % Gamma_D: lower semicircle including (1,0) and (-1,0).
    % Gamma_N: upper semicircle excluding those interface points.
    % The tolerance handles sin(pi), which is not exactly zero numerically.
    tol = 1e-12;
    isDirichlet = isBoundary & (P(:,2) <= tol);
    isNeumann = isBoundary & (P(:,2) > tol);

    [V,A,triangleAreas,neumannEdges] = compute_vertex_weights( ...
        P,T,boundaryEdges,isNeumann);

    E = [T(:,[1 2]); T(:,[2 3]); T(:,[3 1])];
    adjacency = sparse([E(:,1);E(:,2)], [E(:,2);E(:,1)], 1, n, n);
    adjacency = spones(adjacency);

    mesh.points = P;
    mesh.triangles = T;
    mesh.delaunay = DT;
    mesh.boundaryEdges = boundaryEdges;
    mesh.neumannEdges = neumannEdges;
    mesh.isBoundary = isBoundary;
    mesh.isDirichlet = isDirichlet;
    mesh.isNeumann = isNeumann;
    mesh.vertexWeights = V;
    mesh.boundaryWeights = A;
    mesh.triangleAreas = triangleAreas;
    mesh.adjacency = adjacency;
    mesh.area = sum(triangleAreas);
    mesh.neumannBoundaryLength = sum(A);
end


function [V,A,triangleAreas,neumannEdges] = compute_vertex_weights( ...
    P,T,boundaryEdges,isNeumann)
%COMPUTE_VERTEX_WEIGHTS Algorithm 1.

    n = size(P,1);
    e1 = P(T(:,2),:) - P(T(:,1),:);
    e2 = P(T(:,3),:) - P(T(:,1),:);
    triangleAreas = 0.5*abs(e1(:,1).*e2(:,2) - e1(:,2).*e2(:,1));
    V = accumarray(T(:), repmat(triangleAreas/3,3,1), [n,1], @sum, 0);

    % Retain the boundary quadrature rule used in the supplied code:
    % an upper edge shares its full length among its Neumann endpoints.
    midpointY = 0.5*(P(boundaryEdges(:,1),2) + P(boundaryEdges(:,2),2));
    neumannMask = midpointY > 1e-10;
    neumannEdges = boundaryEdges(neumannMask,:);
    A = zeros(n,1);

    for e = 1:size(neumannEdges,1)
        i = neumannEdges(e,1);
        j = neumannEdges(e,2);
        edgeLength = norm(P(i,:) - P(j,:));
        edgeVertices = [i,j];
        edgeNeumannVertices = edgeVertices(isNeumann(edgeVertices));
        if ~isempty(edgeNeumannVertices)
            share = edgeLength/numel(edgeNeumannVertices);
            A(edgeNeumannVertices) = A(edgeNeumannVertices) + share;
        end
    end
end


function sol = solve_discrete_problem(mesh,lambda,mu)
%SOLVE_DISCRETE_PROBLEM Assemble and solve Algorithm 2.

    P = mesh.points;
    n = size(P,1);
    K = assemble_global_matrix(P,mesh.triangles,mesh.vertexWeights, ...
                               mesh.adjacency,lambda,mu);

    x = P(:,1);
    y = P(:,2);
    f = [-16*ones(n,1), zeros(n,1)];
    t = [12*x.^2 + 4*y.^2, 16*x.*y];
    ubar = [ones(n,1), 2*x.*y];

    bNode = bsxfun(@times,mesh.vertexWeights,f) + ...
            bsxfun(@times,mesh.boundaryWeights,t);
    b = reshape(bNode.',[],1);

    D = find(mesh.isDirichlet);
    F = find(~mesh.isDirichlet);
    dofD = nodes_to_dofs(D);
    dofF = nodes_to_dofs(F);

    u = zeros(2*n,1);
    u(dofD) = reshape(ubar(D,:).',[],1);
    KFF = K(dofF,dofF);
    KFD = K(dofF,dofD);
    rhs = b(dofF) - KFD*u(dofD);

    % Always use the dense Moore--Penrose pseudoinverse.
    uFree = pinv(full(KFF))*rhs;
    u(dofF) = uFree;
    U = reshape(u,2,n).';

    Uexact = [x.^2+y.^2, 2*x.*y];
    errorNodes = U - Uexact;
    errorMagnitude = sqrt(sum(errorNodes.^2,2));
    errorMagnitude_mean=mean(errorMagnitude);
    weightedL2 = sqrt(sum(mesh.vertexWeights.*sum(errorNodes.^2,2)));
    weightedL2Exact = sqrt(sum(mesh.vertexWeights.*sum(Uexact.^2,2)));

    % Reconstruct the displacement gradient with exactly the same local
    % coefficients used to assemble K. Rows index displacement components,
    % columns index reference coordinates: G(i,j,a) = partial_j u_i(P_a).
    G = reconstruct_nodal_gradients(P,U,mesh.adjacency);
    Gexact = zeros(2,2,n);
    Gexact(1,1,:) = reshape(2*x,1,1,n);
    Gexact(1,2,:) = reshape(2*y,1,1,n);
    Gexact(2,1,:) = reshape(2*y,1,1,n);
    Gexact(2,2,:) = reshape(2*x,1,1,n);
    gradientError = G - Gexact;
    e11 = reshape(gradientError(1,1,:),n,1);
    e12 = reshape(gradientError(1,2,:),n,1);
    e21 = reshape(gradientError(2,1,:),n,1);
    e22 = reshape(gradientError(2,2,:),n,1);
    diagonalError = sqrt(e11.^2 + e22.^2);
    offDiagonalError = sqrt(e12.^2 + e21.^2);
    fullGradientError = sqrt(diagonalError.^2 + offDiagonalError.^2);

    % This check compares the computed Dirichlet values with the actual
    % prescribed data, rather than with an independently sampled field.
    dirichletDataError = U(D,:) - ubar(D,:);

    sol.K = K;
    sol.load = b;
    sol.bodyForce = f;
    sol.traction = t;
    sol.dirichletData = ubar;
    sol.nodalDisplacement = U;
    sol.exactDisplacement = Uexact;
    sol.analyticalDeformation = P + Uexact;
    sol.discreteDeformation = P + U;
    sol.nodalError = errorNodes;
    sol.nodalErrorMagnitude = errorMagnitude;
    sol.nodalErrorMagnitudeMean = errorMagnitude_mean;
    sol.maxNodalError = max(errorMagnitude);
    sol.weightedL2Error = weightedL2;
    sol.relativeWeightedL2Error = weightedL2/(1+weightedL2Exact);
    sol.dirichletError = max(errorMagnitude(D));
    sol.freeResidual = norm(KFF*uFree-rhs)/(1+norm(rhs));
    sol.dirichletVertices = D;
    sol.freeVertices = F;
    sol.dirichletConstraintError = max(sqrt(sum(dirichletDataError.^2,2)));
    sol.prescribedDirichletDeformation = P(D,:) + ubar(D,:);
    sol.discreteGradient = G;
    sol.exactGradient = Gexact;
    sol.gradientError = gradientError;
    sol.nodalGradientError = fullGradientError;
    sol.nodalDiagonalGradientError = diagonalError;
    sol.nodalOffDiagonalGradientError = offDiagonalError;
    sol.meanGradientError = mean(fullGradientError);
    sol.meanDiagonalGradientError = mean(diagonalError);
    sol.meanOffDiagonalGradientError = mean(offDiagonalError);
    sol.gradientComponentMeanAbsoluteError = ...
        [mean(abs(e11)), mean(abs(e12)); ...
         mean(abs(e21)), mean(abs(e22))];
end


function K = assemble_global_matrix(P,T,V,adjacency,lambda,mu)
%ASSEMBLE_GLOBAL_MATRIX Assemble Equation (22) in sparse block form.

    n = size(P,1);
    I2 = eye(2);
    rowParts = cell(n,1);
    colParts = cell(n,1);
    valParts = cell(n,1);

    for j = 1:n
        [localVertices,localP] = local_reconstruction(P,adjacency,j);
        nlocal = numel(localVertices);

        numberOfEntries = 4*nlocal*nlocal;
        I = zeros(numberOfEntries,1);
        J = zeros(numberOfEntries,1);
        S = zeros(numberOfEntries,1);
        counter = 0;

        for a = 1:nlocal
            p_i = localP(a,:);
            dofI = [2*localVertices(a)-1, 2*localVertices(a)];
            for q = 1:nlocal
                p_q = localP(q,:);
                dofQ = [2*localVertices(q)-1, 2*localVertices(q)];
                block = mu*(p_i*p_q.')*I2 + ...
                        mu*(p_q.')*p_i + lambda*(p_i.')*p_q;
                block = V(j)*block;
                for r = 1:2
                    for c = 1:2
                        counter = counter + 1;
                        I(counter) = dofI(r);
                        J(counter) = dofQ(c);
                        S(counter) = block(r,c);
                    end
                end
            end
        end

        rowParts{j} = I;
        colParts{j} = J;
        valParts{j} = S;
    end

    I = vertcat(rowParts{:});
    J = vertcat(colParts{:});
    S = vertcat(valParts{:});
    K = sparse(I,J,S,2*n,2*n);
    K = 0.5*(K + K.');
end


function [localVertices,localP] = local_reconstruction(P,adjacency,j)
%LOCAL_RECONSTRUCTION Shared least-squares coefficients for assembly and G.

    neighbours = find(adjacency(:,j));
    differences = P(neighbours,:) - repmat(P(j,:),numel(neighbours),1);
    C = differences.'*differences;
    R = differences*pinv(C);
    localVertices = [neighbours; j];
    localP = [R; -sum(R,1)];
end


function G = reconstruct_nodal_gradients(P,U,adjacency)
%RECONSTRUCT_NODAL_GRADIENTS Return a 2-by-2-by-n displacement gradient.
% G(:,:,j) = DeltaU_j.' * X_j * pinv(X_j.' * X_j), equivalently
%           U(localVertices,:).' * localP.

    n = size(P,1);
    G = zeros(2,2,n);
    for j = 1:n
        [localVertices,localP] = local_reconstruction(P,adjacency,j);
        G(:,:,j) = U(localVertices,:).'*localP;
    end
end


function dof = nodes_to_dofs(nodes)
%NODES_TO_DOFS Convert vertex indices to interleaved [x_i,y_i] DOFs.

    nodes = nodes(:);
    dof = reshape([2*nodes-1, 2*nodes].',[],1);
end


function print_summary(name,mesh,sol)
    fprintf('%s discretization:\n',name);
    fprintf('  vertices             : %d\n',size(mesh.points,1));
    fprintf('  triangles            : %d\n',size(mesh.triangles,1));
    fprintf('  Dirichlet vertices   : %d\n',sum(mesh.isDirichlet));
    fprintf('  Neumann vertices     : %d\n',sum(mesh.isNeumann));
    fprintf('  L2 Error Mesh   : %.6e\n',sol.relativeWeightedL2Error);
    fprintf('  Mean error over vertices : %.6e\n',sol.nodalErrorMagnitudeMean);
    fprintf('  mean gradient error   : %.6e\n',sol.meanGradientError);
    fprintf('  mean diagonal error   : %.6e\n',sol.meanDiagonalGradientError);
    fprintf('  mean off-diagonal err.: %.6e\n',sol.meanOffDiagonalGradientError);
    M = sol.gradientComponentMeanAbsoluteError;
    fprintf('  Entrywise mean absolute gradient errors (row i, column j):\n');
    fprintf('    [%.6e  %.6e]\n',M(1,1),M(1,2));
    fprintf('    [%.6e  %.6e]\n\n',M(2,1),M(2,2));
end


function plot_meshes(meshRand,meshEq)
%PLOT_MESHES Visualize the two point sets and their triangulations.

    figure('Color','w','Name','Delaunay discretizations');

    subplot(1,2,1);
    triplot(meshRand.triangles,meshRand.points(:,1),meshRand.points(:,2), ...
            'Color',[0.75 0.75 0.75]);
    hold on;
    plot(meshRand.points(meshRand.isDirichlet,1), ...
         meshRand.points(meshRand.isDirichlet,2),'r.','MarkerSize',10);
    plot(meshRand.points(meshRand.isNeumann,1), ...
         meshRand.points(meshRand.isNeumann,2),'b.','MarkerSize',10);
    axis equal tight;
    title('Random point distribution');
    xlabel('x'); ylabel('y');
   legend('Delaunay edges', '\Gamma_D', '\Gamma_N', 'Location', 'southoutside');

    subplot(1,2,2);
    triplot(meshEq.triangles,meshEq.points(:,1),meshEq.points(:,2), ...
            'Color',[0.75 0.75 0.75]);
    hold on;
    plot(meshEq.points(meshEq.isDirichlet,1), ...
         meshEq.points(meshEq.isDirichlet,2),'r.','MarkerSize',10);
    plot(meshEq.points(meshEq.isNeumann,1), ...
         meshEq.points(meshEq.isNeumann,2),'b.','MarkerSize',10);
    axis equal tight;
    title('Quasi-uniform point distribution');
    xlabel('x'); ylabel('y');
    legend('Delaunay edges','\Gamma_D','\Gamma_N','Location','southoutside');
end


function plot_solutions(meshRand,solRand,meshEq,solEq)
%PLOT_SOLUTIONS Plot the first displacement component and its nodal error.

    % Create a wider figure window with a custom size [left, bottom, width, height]
    figure('Color', 'w', 'Name', 'Discrete error', 'Position', [100, 100, 1200, 550]);

    % Plot random mesh error on the left subplot
    subplot(1, 2, 1);
    trisurf(meshRand.triangles, meshRand.points(:,1), meshRand.points(:,2), ...
            solRand.nodalErrorMagnitude, 'EdgeColor', 'none');
    view(2); axis equal tight; colorbar;
    title(sprintf('Random error, max = %.2e', solRand.maxNodalError));
    xlabel('x'); ylabel('y');

    % Plot quasi-uniform mesh error on the right subplot
    subplot(1, 2, 2);
    trisurf(meshEq.triangles, meshEq.points(:,1), meshEq.points(:,2), ...
            solEq.nodalErrorMagnitude, 'EdgeColor', 'none');
    view(2); axis equal tight; colorbar;
    title(sprintf('Quasi-uniform error, max = %.2e', solEq.maxNodalError));
    xlabel('x'); ylabel('y');
end


function plot_deformations(meshRand,solRand,meshEq,solEq)
%PLOT_DEFORMATIONS Group six point-cloud views in one large 2-by-3 figure.
% Rows are random/quasi-uniform. Columns are discrete deformation,
% original disk, and analytical deformation alone. All views retain the
% same reference-node colors and use Phi = X + U with displacement scale 1.
% The color key uses the original angle theta_0 and radius r_0, never the
% deformed coordinates. Six configuration panels remain in a 2-by-3 layout.

    meshes = {meshRand,meshEq};
    solutions = {solRand,solEq};
    names = {'Random','Quasi-uniform'};

    % Use common limits for the two analytical-only panels.
    allAnalytical = [solRand.analyticalDeformation;solEq.analyticalDeformation];
    lower = min(allAnalytical,[],1);
    upper = max(allAnalytical,[],1);
    span = upper - lower;
    padding = 0.05*span;
    padding(padding == 0) = 0.05*max([span,1]);
    analyticalXLimits = [lower(1)-padding(1),upper(1)+padding(1)];
    analyticalYLimits = [lower(2)-padding(2),upper(2)+padding(2)];

    create_large_figure('Six views: original disks and deformations');
    for k = 1:2
        mesh = meshes{k};
        sol = solutions{k};
        colors = reference_node_colors(mesh.points);
        firstPanel = 3*(k-1);
        % Leave a bottom strip for the two-dimensional color key.
        panelBottom = 0.61 - 0.37*(k-1);

        subplot(2,3,firstPanel+1);
        set(gca,'Position',[0.055,panelBottom,0.27,0.28]);
        discrete = sol.discreteDeformation;
        lowerDiscrete = min(discrete,[],1);
        upperDiscrete = max(discrete,[],1);
        spanDiscrete = upperDiscrete - lowerDiscrete;
        paddingDiscrete = 0.05*spanDiscrete;
        paddingDiscrete(paddingDiscrete == 0) = 0.05*max([spanDiscrete,1]);
        draw_colored_point_cloud(discrete,colors, ...
            [lowerDiscrete(1)-paddingDiscrete(1),upperDiscrete(1)+paddingDiscrete(1)], ...
            [lowerDiscrete(2)-paddingDiscrete(2),upperDiscrete(2)+paddingDiscrete(2)], ...
            mesh.isBoundary);
        title({names{k},'Discrete deformation'},'FontSize',12);

        subplot(2,3,firstPanel+2);
        set(gca,'Position',[0.370,panelBottom,0.27,0.28]);
        draw_colored_point_cloud(mesh.points,colors, ...
            [-1.05,1.05],[-1.05,1.05],mesh.isBoundary);
        title({names{k},'Original disk'},'FontSize',12);

        subplot(2,3,firstPanel+3);
        set(gca,'Position',[0.685,panelBottom,0.27,0.28]);
        draw_colored_point_cloud(sol.analyticalDeformation,colors, ...
            analyticalXLimits,analyticalYLimits,mesh.isBoundary);
        title({names{k},'Analytical deformation'},'FontSize',12);
    end
    draw_reference_color_key;
    annotation('textbox',[0.055 0.035 0.24 0.11], ...
        'String',{'Hue: original angle.', ...
                  'Pale to intense/darker: original radius.', ...
                  'Each node keeps its color in all three views.'}, ...
        'HorizontalAlignment','left','EdgeColor','none','FontSize',10);
annotation('textbox',[0.745 0.035 0.22 0.11], ...
    'String','Outlined markers: boundary nodes.', ...
    'HorizontalAlignment','left','EdgeColor','none','FontSize',10);
end


function draw_colored_point_cloud(points,colors,xLimits,yLimits,isBoundary)
%DRAW_COLORED_POINT_CLOUD Plot a single configuration with reference colors.
% Draw interiors first so that boundary nodes remain visible in both meshes.

    markerArea = max(8,min(32,32000/size(points,1)));
    interior = ~isBoundary;
    scatter(points(interior,1),points(interior,2),markerArea, ...
        colors(interior,:),'o','filled', ...
        'MarkerEdgeColor','none');
    hold on;
    scatter(points(isBoundary,1),points(isBoundary,2),1.5*markerArea, ...
        colors(isBoundary,:),'o','filled', ...
        'MarkerEdgeColor',[0.20 0.20 0.20],'LineWidth',0.5);
    hold off;
    axis equal;
    xlim(xLimits);
    ylim(yLimits);
    set(gca,'FontSize',12);
   % xlabel('x'); ylabel('y');
    box on;
end


function colors = reference_node_colors(P)
%REFERENCE_NODE_COLORS Assign colors using reference angle and radius.
% The rule is identical for both meshes, including their shared boundary.
% The original position is represented by its polar coordinates (theta,r).
% Increasing r gives stronger saturation and lower brightness, making
% points on the same ray distinguishable. Never recolor after deformation.

    angle = mod(atan2(P(:,2),P(:,1)),2*pi)/(2*pi);
    radius = min(1,sqrt(sum(P.^2,2)));
    saturation = 0.20 + 0.75*radius;
    brightness = 0.98 - 0.28*radius;
    % The angle is undefined at the center: use a neutral gray there.
    atCenter = radius <= 1e-12;
    saturation(atCenter) = 0;
    brightness(atCenter) = 0.85;
    colors = hsv2rgb([angle,saturation,brightness]);
end


function draw_reference_color_key
%DRAW_REFERENCE_COLOR_KEY Display the exact bivariate reference color rule.
% Truecolor RGB is used, so the key does not affect other figure colormaps.

    angleDegrees = linspace(0,360,361);
    radii = linspace(0,1,101);
    [Theta,Radius] = meshgrid(angleDegrees*pi/180,radii);
    referencePoints = [Radius(:).*cos(Theta(:)), ...
                       Radius(:).*sin(Theta(:))];
    rgb = reference_node_colors(referencePoints);
    rgbImage = reshape(rgb,[numel(radii),numel(angleDegrees),3]);

    keyAxes = axes('Position',[0.365 0.060 0.345 0.075]);
    image(keyAxes,angleDegrees,radii,rgbImage);
    set(keyAxes,'YDir','normal','XTick',0:90:360, ...
        'YTick',[0,0.5,1],'FontSize',9);
    xlim(keyAxes,[0,360]);
    ylim(keyAxes,[0,1]);
    xlabel(keyAxes,'Original angle \theta_0 (degrees)');
    ylabel(keyAxes,'Radius r_0');
    title(keyAxes,'Color key: original position','FontSize',10);
    box(keyAxes,'on');
end




function plot_gradient_errors(meshRand,solRand,meshEq,solEq)
%PLOT_GRADIENT_ERRORS Group all six gradient-error maps in a 3-by-2 figure.
% Each row retains common color limits for the two discretizations.

    fieldsRand = {solRand.nodalGradientError, ...
                  solRand.nodalDiagonalGradientError, ...
                  solRand.nodalOffDiagonalGradientError};
    fieldsEq = {solEq.nodalGradientError, ...
                solEq.nodalDiagonalGradientError, ...
                solEq.nodalOffDiagonalGradientError};
    labels = {'Full gradient error','Diagonal error','Off-diagonal error'};

    create_large_figure('Nodal displacement-gradient errors');
    for r = 1:3
        maxError = max([fieldsRand{r};fieldsEq{r}]);
        if maxError <= 0
            maxError = 1;
        end

        subplot(3,2,2*r-1);
        draw_gradient_error(meshRand,fieldsRand{r},maxError);
        title({['Random: ' labels{r}], ...
            sprintf('Mean = %.2e',mean(fieldsRand{r}))},'FontSize',12);

        subplot(3,2,2*r);
        draw_gradient_error(meshEq,fieldsEq{r},maxError);
        title({['Quasi-uniform: ' labels{r}], ...
            sprintf('Mean = %.2e',mean(fieldsEq{r}))},'FontSize',12);
    end
end


function draw_gradient_error(mesh,errorValues,maxError)
%DRAW_GRADIENT_ERROR Display a nodal error field on the reference mesh.

    trisurf(mesh.triangles,mesh.points(:,1),mesh.points(:,2),errorValues, ...
        'FaceColor','interp','EdgeColor','none');
    view(2);
    axis equal tight;
    caxis([0,maxError]);
    colorbar;
    set(gca,'FontSize',11);
    xlabel('x'); ylabel('y');
end


function create_large_figure(name)
%CREATE_LARGE_FIGURE Use most of the screen without fixed pixel dimensions.
% Normalized placement and ordinary windows are compatible with R2017b.

    figure('Color','w','Name',name,'NumberTitle','off', ...
        'WindowStyle','normal','Units','normalized', ...
        'OuterPosition',[0.01 0.04 0.98 0.92]);
end
