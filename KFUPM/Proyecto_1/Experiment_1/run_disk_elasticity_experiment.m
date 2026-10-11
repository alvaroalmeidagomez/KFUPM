function results = run_disk_elasticity_experiment(Nlevels)
% =========================================================================
% Developed at King Fahd University of Petroleum and Minerals (KFUPM)
% Dhahran, Saudi Arabia -- IRC for Industrial Nuclear Energy
% =========================================================================
% Vertex-based stationary elasticity: random versus regular interior meshes.
% Compatible with MATLAB R2017b. Save under this exact function filename.
%
% Run interactively:
%   results = run_disk_elasticity_experiment;
% Or supply target node counts:
%   results = run_disk_elasticity_experiment([100 400 1600]);
% A scalar N runs up to three levels: max(30,round(N/4)), max(30,round(N/2)), N.
% Repeated counts are removed. At N=30 only one level is available (p = NaN).
% Every free system is solved with pinv(full(KFF)); plots are always displayed.
% Dense pseudoinverses become expensive as node counts increase.
%
% The regular family uses a Cartesian lattice in the interior and a transition
% to the same equally spaced circular boundary as the random family. The actual
% regular node count may differ slightly from its target to retain the lattice.
% Complete interior squares have identical bottom-left to top-right diagonals,
% imposed as constraints in the Delaunay triangulation.
% Boundary count is approximately 10 percent of the TARGET count, rounded up
% to an even number, at least 4, to include both boundary-condition interfaces.
%
% The computational domain is the inscribed boundary polygon. Manufactured
% traction is integrated on its upper straight edges as sigma_exact*n_edge,
% using two-point Gauss quadrature and linear endpoint test values. No endpoint
% weight is reassigned when the other endpoint is Dirichlet. Body-force weights
% and the vertex-based stiffness/reconstruction are retained from the input.
%
% Exact fields (small strain, two-dimensional isotropic constitutive law):
%   u = [x^2+y^2; 2*x*y]
%   epsilon = [2*x, 2*y; 2*y, 2*x]
%   sigma = [4*(lambda+mu)*x, 4*mu*y; 4*mu*y, 4*(lambda+mu)*x]
%   f = [-4*lambda-8*mu; 0].  Default lambda=2, mu=1.
% No plane-stress parameter reduction is applied.
%
% Nodal displacement errors use vector Euclidean norms. Tensor errors use the
% full Frobenius norm, including both off-diagonal entries. Weighted L2 denotes
% sqrt(sum_a V_a*error_a^2), a nodal quadrature norm, not an exact integral.
% True relative errors divide by the corresponding nonzero exact field norm.
% h is the maximum Delaunay edge length on each actual triangulation.
% p_obs = log(E_coarse/E_fine)/log(h_coarse/h_fine), computed separately for each
% family, for both weighted L2 and maximum errors. Invalid orders are NaN.
% Random rates use one seeded realization per level; they can fluctuate.
% Observed orders may be negative; no positive convergence rate is assumed.
%
% Free equilibrium residuals K*u-b and K*u_exact-b are stored separately.
% Relative residuals use norm(b_F-K_FD*u_D), without an added constant.
% Raw nodal force residuals are not pointwise strong-form PDE residuals.
% Displacement-gradient reconstruction errors on exact samples are also stored.
% Gradient comparison uses E_G = G_h - grad(u_exact), with G_ij = d_j u_i.
% Diagonal errors group E11,E22; off-diagonal errors group E12,E21. The console
% reports each component and each group at EVERY level, for both mesh families.
% Component tables contain mean absolute, maximum absolute, weighted L2 and
% relative L2 errors; group tables use Euclidean norms of the selected entries.
% Signed nodal components are stored in solution.gradientComponentErrors,
% in the order G11,G22,G12,G21. Summary tables are retained with each solution
% and in results.convergence.<family>.gradientErrorTable.
%
% Figures use the finest level. Full fields and diagnostics from EVERY level
% are retained in results.levels. Tables are in results.convergence.
% The function saves results to Data_files.mat in the current MATLAB folder.


    clc;
    close all;
    
disp('=========================================================================');
disp('Developed at King Fahd University of Petroleum and Minerals (KFUPM)');
disp('Dhahran, Saudi Arabia');
disp('IRC for Industrial Nuclear Energy');
disp('=========================================================================');

    if nargin < 1
        clc;
        disp('KFUPM -- Vertex-based elasticity: random and regular meshes');
        while true
            textN = input(['Enter target counts, e.g. [100 400 1600], ' ...
                           'or a finest count N >= 30: '], 's');
            % Numeric parsing only; never evaluate an input expression.
            cleaned = regexprep(textN,'[\[\],;]',' ');
            if ~isempty(regexp(cleaned,'^[0-9eE+\-.\s]+$','once'))
                tokens = regexp(strtrim(cleaned),'\s+','split');
                Nlevels = str2double(tokens);
            else
                Nlevels = [];
            end
            if valid_counts(Nlevels)
                break;
            end
            fprintf('Use integer counts >= 30, separated by spaces or commas.\n');
        end
    elseif ~valid_counts(Nlevels)
        error('Node counts must be a nonempty real vector of integers >= 30.');
    end

    Nlevels = double(Nlevels(:).');
    if numel(Nlevels) == 1
        finest = Nlevels;
        Nlevels = [max(30,round(finest/4)), max(30,round(finest/2)), finest];
    end
    Nlevels = unique(Nlevels); % Ascending target counts.
    lambda = 2;
    mu = 1;
    numberOfLevels = numel(Nlevels);
    results.parameters.targetNodeCounts = Nlevels;
    results.parameters.lambda = lambda;
    results.parameters.mu = mu;
    results.parameters.makePlots = true;
    results.parameters.useMoorePenrose = true;
    results.parameters.deformationScale = 1;
    results.parameters.hDefinition = 'Maximum unique Delaunay edge length';
    results.parameters.boundaryLoadRule = 'Polygon normals; two-point Gauss';
    results.parameters.randomSeed = 1;
    results.parameters.regularMesh = 'Cartesian interior, fixed diagonals, polygon boundary';
    rng(1,'twister');

    for level = 1:numberOfLevels
        Ntarget = Nlevels(level);
        Nboundary = max(4,2*ceil(round(0.10*Ntarget)/2));
        theta = (0:Nboundary-1)'*(2*pi/Nboundary);
        boundaryPoints = [cos(theta),sin(theta)];
        P_rand = make_random_point_set(Ntarget,boundaryPoints);
        [P_regular,regularInfo] = make_regular_point_set(Ntarget,boundaryPoints);
        meshRand = prepare_mesh(P_rand);
        meshRegular = prepare_mesh(P_regular,regularInfo.constraints);
        meshRand.targetNodeCount = Ntarget;
        meshRegular.targetNodeCount = Ntarget;
        meshRegular.regularGridSpacing = regularInfo.gridSpacing;
        meshRegular.regularGridInfo = regularInfo;

        fprintf('\nLEVEL %d/%d: target N = %d; shared boundary = %d\n', ...
            level,numberOfLevels,Ntarget,Nboundary);
        solRand = solve_discrete_problem(meshRand,lambda,mu);
        solRegular = solve_discrete_problem(meshRegular,lambda,mu);
        print_summary('Random',meshRand,solRand);
        print_summary('Regular',meshRegular,solRegular);

        results.levels(level).targetNodeCount = Ntarget;
        results.levels(level).boundaryPoints = boundaryPoints;
        results.levels(level).random.mesh = meshRand;
        results.levels(level).random.solution = solRand;
        results.levels(level).regular.mesh = meshRegular;
        results.levels(level).regular.solution = solRegular;
    end

    results.random = results.levels(end).random;
    results.regular = results.levels(end).regular;
    results.parameters.N = Nlevels(end); % Finest TARGET; see actual mesh counts.
    results.parameters.Nboundary = size(results.levels(end).boundaryPoints,1);
    results.parameters.boundaryPoints = results.levels(end).boundaryPoints;
    results.parameters.boundaryAngleStep = 2*pi/results.parameters.Nboundary;
    results.convergence.random = build_convergence(results.levels,'random');
    results.convergence.regular = build_convergence(results.levels,'regular');
    print_convergence('Random',results.convergence.random);
    print_convergence('Regular',results.convergence.regular);
    if numberOfLevels < 2
        fprintf('Only one level: observed orders are undefined (NaN).\n');
    end

    meshRand = results.random.mesh;
    meshRegular = results.regular.mesh;
    solRand = results.random.solution;
    solRegular = results.regular.solution;
    plot_meshes(meshRand,meshRegular);
    plot_field_comparisons(meshRand,solRand,meshRegular,solRegular);
    plot_deformations(meshRand,solRand,meshRegular,solRegular);
    plot_gradient_errors(meshRand,solRand,meshRegular,solRegular);
    plot_free_residuals(meshRand,solRand,meshRegular,solRegular);
    plot_convergence(results.convergence);
    save('Data_files.mat','results','-v7.3');
    fprintf('\nAll levels and diagnostic tables saved to Data_files.mat.\n');
end


function valid = valid_counts(counts)
    valid = isnumeric(counts) && isreal(counts) && isvector(counts) && ...
        ~isempty(counts) && all(isfinite(counts(:))) && ...
        all(counts(:) == fix(counts(:))) && all(counts(:) >= 30);
end


function P = make_random_point_set(N,boundaryPoints)
    Nboundary = size(boundaryPoints,1);
    Ninterior = N - Nboundary;
    rmax = cos(pi/Nboundary)*(1-1e-12);
    interior = zeros(Ninterior,2);
    nfilled = 0;
    while nfilled < Ninterior
        batchSize = max(100,2*(Ninterior-nfilled));
        candidate = 2*rand(batchSize,2)-1;
        candidate = candidate(sum(candidate.^2,2) < rmax^2,:);
        nnew = min(size(candidate,1),Ninterior-nfilled);
        interior(nfilled+(1:nnew),:) = candidate(1:nnew,:);
        nfilled = nfilled+nnew;
    end
    P = [boundaryPoints;interior];
end


function [P,meta] = make_regular_point_set(N,boundaryPoints)
%MAKE_REGULAR_POINT_SET Cartesian lattice with the supplied polygon boundary.
% Compatible with MATLAB R2017b. The boundary points are retained exactly,
% in their supplied order, and appear first in P. Interior points have the
% form (i*h,j*h), for integer i,j, with one common grid spacing h.
%
% N is a TARGET total node count. Keeping a complete regular lattice can
% prevent attaining exactly N; meta.Nactual reports the actual count.
% A clearance of 0.25*h from each polygon edge avoids nearly coincident
% boundary and lattice nodes. Consequently, the Delaunay mesh has a regular
% Cartesian interior and an unavoidable transition next to the circle.
% Use delaunayTriangulation(P,meta.constraints) to split every complete
% Cartesian square along the same bottom-left to top-right diagonal.

    if ~isscalar(N) || ~isfinite(N) || N ~= floor(N)
        error('N must be a finite integer.');
    end
    if size(boundaryPoints,2) ~= 2 || size(boundaryPoints,1) < 3 || ...
            any(~isfinite(boundaryPoints(:)))
        error('boundaryPoints must contain at least three finite 2D points.');
    end
    nb = size(boundaryPoints,1);
    targetInterior = N - nb;
    if targetInterior < 1
        error('N must exceed the number of boundary points.');
    end

    following = boundaryPoints([2:nb,1],:);
    edges = following - boundaryPoints;
    edgeLengths = sqrt(sum(edges.^2,2));
    if any(edgeLengths == 0)
        error('Consecutive boundary points must be distinct.');
    end
    signedArea = 0.5*sum(boundaryPoints(:,1).*following(:,2) - ...
                         boundaryPoints(:,2).*following(:,1));
    if signedArea == 0
        error('The boundary polygon must have positive area.');
    end
    orientation = sign(signedArea);
    nextEdges = edges([2:nb,1],:);
    turns = orientation*(edges(:,1).*nextEdges(:,2) - ...
                         edges(:,2).*nextEdges(:,1));
    if any(turns < -1e-12)
        error('The boundary points must describe a convex polygon in order.');
    end
    [originInside,originOn] = inpolygon(0,0, ...
        boundaryPoints(:,1),boundaryPoints(:,2));
    if ~originInside || originOn
        error('The origin must lie strictly inside the boundary polygon.');
    end

    clearanceFraction = 0.25;
    h = sqrt(abs(signedArea)/targetInterior);
    interior = regular_lattice_in_polygon( ...
        h,boundaryPoints,edges,edgeLengths,orientation,clearanceFraction);
    count = size(interior,1);
    bestH = h;
    bestInterior = interior;
    bestDifference = abs(count - targetInterior);

    % The count is nonincreasing in h: the polygon contains the origin and
    % both each lattice point and its required clearance scale linearly.
    low = h;
    high = h;
    if count > targetInterior
        while count > targetInterior
            high = 1.5*high;
            interior = regular_lattice_in_polygon(high,boundaryPoints, ...
                edges,edgeLengths,orientation,clearanceFraction);
            count = size(interior,1);
        end
        h = high;
    elseif count < targetInterior
        while count < targetInterior
            low = low/1.5;
            interior = regular_lattice_in_polygon(low,boundaryPoints, ...
                edges,edgeLengths,orientation,clearanceFraction);
            count = size(interior,1);
        end
        h = low;
    end
    if abs(count-targetInterior) < bestDifference
        bestInterior = interior;
        bestDifference = abs(count-targetInterior);
        bestH = h;
    end

    for iteration = 1:28
        if bestDifference == 0 || high == low
            break;
        end
        h = 0.5*(low + high);
        interior = regular_lattice_in_polygon(h,boundaryPoints,edges, ...
            edgeLengths,orientation,clearanceFraction);
        count = size(interior,1);
        difference = abs(count - targetInterior);
        if difference < bestDifference
            bestH = h;
            bestInterior = interior;
            bestDifference = difference;
        end
        if count >= targetInterior
            low = h;
        else
            high = h;
        end
    end

    P = [boundaryPoints; bestInterior];
    meta.Ntarget = N;
    meta.Nactual = size(P,1);
    meta.Nboundary = nb;
    meta.Ninterior = size(bestInterior,1);
    meta.gridSpacing = bestH;
    meta.gridOrigin = [0,0];
    meta.boundaryClearance = clearanceFraction*bestH;
    meta.clearanceFraction = clearanceFraction;
    meta.polygonArea = abs(signedArea);
    [meta.constraints,meta.regularTriangles] = ...
        regular_grid_constraints(P,nb,bestH);
    meta.NregularSquares = size(meta.regularTriangles,1)/2;
end

function [constraints,regularTriangles] = ...
    regular_grid_constraints(P,Nboundary,spacing)
%REGULAR_GRID_CONSTRAINTS Repeating triangulation of complete lattice cells.
% Every complete square has its four sides and the bottom-left to top-right
% diagonal constrained. Boundary polygon edges are constrained as well.
% regularTriangles records the two prescribed faces of each full square;
% faces in the boundary transition are supplied by delaunayTriangulation.

    boundaryEdges = [(1:Nboundary)', [2:Nboundary,1]'];
    interiorIndices = (Nboundary+1:size(P,1))';
    if isempty(interiorIndices)
        constraints = boundaryEdges;
        regularTriangles = zeros(0,3);
        return;
    end
    integerCoordinates = round(P(interiorIndices,:)/spacing);
    minX = min(integerCoordinates(:,1));
    minY = min(integerCoordinates(:,2));
    ncolumns = max(integerCoordinates(:,1)) - minX + 1;
    nrows = max(integerCoordinates(:,2)) - minY + 1;
    grid = zeros(nrows,ncolumns);
    rows = integerCoordinates(:,2) - minY + 1;
    columns = integerCoordinates(:,1) - minX + 1;
    grid(sub2ind([nrows,ncolumns],rows,columns)) = interiorIndices;

    bottomLeft = grid(1:end-1,1:end-1);
    bottomRight = grid(1:end-1,2:end);
    topLeft = grid(2:end,1:end-1);
    topRight = grid(2:end,2:end);
    fullSquare = bottomLeft > 0 & bottomRight > 0 & ...
                 topLeft > 0 & topRight > 0;
    a = bottomLeft(fullSquare);
    b = bottomRight(fullSquare);
    c = topLeft(fullSquare);
    d = topRight(fullSquare);
    a = a(:); b = b(:); c = c(:); d = d(:);
    squareEdges = [a,b; a,c; b,d; c,d; a,d];
    constraints = unique(sort([boundaryEdges; squareEdges],2),'rows');
    regularTriangles = [a,b,d; a,d,c];
end

function interior = regular_lattice_in_polygon( ...
    h,boundaryPoints,edges,edgeLengths,orientation,clearanceFraction)
% Retain all lattice points in the inward-offset convex polygon.
    x = (ceil(min(boundaryPoints(:,1))/h): ...
         floor(max(boundaryPoints(:,1))/h))*h;
    y = (ceil(min(boundaryPoints(:,2))/h): ...
         floor(max(boundaryPoints(:,2))/h))*h;
    [X,Y] = meshgrid(x,y);
    interior = [X(:),Y(:)];
    clearance = clearanceFraction*h;
    for edge = 1:size(boundaryPoints,1)
        dx = interior(:,1) - boundaryPoints(edge,1);
        dy = interior(:,2) - boundaryPoints(edge,2);
        signedDistance = orientation*(edges(edge,1)*dy - ...
            edges(edge,2)*dx)/edgeLengths(edge);
        interior = interior(signedDistance >= clearance,:);
    end
end


function mesh = prepare_mesh(P,constraintEdges)
    n = size(P,1);
    if nargin < 2 || isempty(constraintEdges)
        DT = delaunayTriangulation(P);
    else
        DT = delaunayTriangulation(P,constraintEdges);
    end
    T = DT.ConnectivityList;
    boundaryEdges = freeBoundary(DT);
    isBoundary = false(n,1);
    isBoundary(unique(boundaryEdges(:))) = true;
    tol = 1e-12;
    isDirichlet = isBoundary & (P(:,2) <= tol);
    isNeumann = isBoundary & (P(:,2) > tol);
    [V,A,triangleAreas,neumannEdges] = compute_vertex_weights( ...
        P,T,boundaryEdges,isNeumann);
    E = unique(sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2),'rows');
    adjacency = sparse([E(:,1);E(:,2)],[E(:,2);E(:,1)],1,n,n);
    adjacency = spones(adjacency);
    lengths = sqrt(sum((P(E(:,1),:)-P(E(:,2),:)).^2,2));
    mesh.points = P;
    mesh.triangles = T;
    mesh.delaunay = DT;
    mesh.edges = E;
    mesh.edgeLengths = lengths;
    mesh.h = max(lengths);
    mesh.hMean = mean(lengths);
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
% Each triangle contributes one third of its area to each vertex.
% The even boundary count places vertices at both mixed-BC interfaces.
% Thus no boundary edge crosses an interface; select upper edges by endpoints.
    n = size(P,1);
    e1 = P(T(:,2),:)-P(T(:,1),:);
    e2 = P(T(:,3),:)-P(T(:,1),:);
    triangleAreas = 0.5*abs(e1(:,1).*e2(:,2)-e1(:,2).*e2(:,1));
    V = accumarray(T(:),repmat(triangleAreas/3,3,1),[n,1],@sum,0);
    endpointY = P(boundaryEdges(:,1),2);
    otherY = P(boundaryEdges(:,2),2);
    upperEdges = endpointY >= -1e-12 & otherY >= -1e-12 & ...
        (isNeumann(boundaryEdges(:,1)) | isNeumann(boundaryEdges(:,2)));
    neumannEdges = boundaryEdges(upperEdges,:);
    A = zeros(n,1);
    for e = 1:size(neumannEdges,1)
        vertices = neumannEdges(e,:);
        lengthEdge = norm(P(vertices(2),:)-P(vertices(1),:));
        A(vertices) = A(vertices)+lengthEdge/2;
    end
% A includes the half weights at interface nodes. Their test values are zero;
% their load entries disappear when restricting to free DOFs, not by doubling
% their neighbours' weights. Actual vector traction loads are integrated below.
end


function boundaryLoad = integrate_neumann_load(mesh,lambda,mu)
% Exact degree-two endpoint-weighted traction quadrature on polygon edges.
    P = mesh.points;
    boundaryLoad = zeros(size(P));
    quadraturePoints = 0.5+[-1,1]/(2*sqrt(3));
    for e = 1:size(mesh.neumannEdges,1)
        vertices = mesh.neumannEdges(e,:);
        a = P(vertices(1),:);
        b = P(vertices(2),:);
        tangent = b-a;
        edgeLength = norm(tangent);
        normal = [tangent(2);-tangent(1)]/edgeLength;
        if dot(normal,0.5*(a+b).') < 0
            normal = -normal; % Outward for an origin-containing polygon.
        end
        for q = 1:2
            s = quadraturePoints(q);
            point = (1-s)*a+s*b;
            sigma = [4*(lambda+mu)*point(1),4*mu*point(2); ...
                     4*mu*point(2),4*(lambda+mu)*point(1)];
            traction = (sigma*normal).';
            boundaryLoad(vertices(1),:) = boundaryLoad(vertices(1),:) + ...
                (edgeLength/2)*(1-s)*traction;
            boundaryLoad(vertices(2),:) = boundaryLoad(vertices(2),:) + ...
                (edgeLength/2)*s*traction;
        end
    end
end


function sol = solve_discrete_problem(mesh,lambda,mu)
    P = mesh.points;
    n = size(P,1);
    K = assemble_global_matrix(P,mesh.triangles,mesh.vertexWeights, ...
                               mesh.adjacency,lambda,mu);
    x = P(:,1);
    y = P(:,2);
    Uexact = [x.^2+y.^2,2*x.*y];
    f = [(-4*lambda-8*mu)*ones(n,1),zeros(n,1)];
    boundaryLoad = integrate_neumann_load(mesh,lambda,mu);
    bNode = bsxfun(@times,mesh.vertexWeights,f)+boundaryLoad;
    b = reshape(bNode.',[],1);
    D = find(mesh.isDirichlet);
    F = find(~mesh.isDirichlet);
    dofD = nodes_to_dofs(D);
    dofF = nodes_to_dofs(F);
    u = zeros(2*n,1);
    u(dofD) = reshape(Uexact(D,:).',[],1);
    KFF = K(dofF,dofF);
    rhs = b(dofF)-K(dofF,dofD)*u(dofD);
    u(dofF) = pinv(full(KFF))*rhs;
    U = reshape(u,2,n).';
    G = reconstruct_nodal_gradients(P,U,mesh.adjacency);
    Gexact = zeros(2,2,n);
    Gexact(1,1,:) = reshape(2*x,1,1,n);
    Gexact(1,2,:) = reshape(2*y,1,1,n);
    Gexact(2,1,:) = reshape(2*y,1,1,n);
    Gexact(2,2,:) = reshape(2*x,1,1,n);
    strain = 0.5*(G+permute(G,[2,1,3]));
    exactStrain = Gexact;
    stress = 2*mu*strain;
    traceStrain = strain(1,1,:)+strain(2,2,:);
    stress(1,1,:) = stress(1,1,:)+lambda*traceStrain;
    stress(2,2,:) = stress(2,2,:)+lambda*traceStrain;
    exactStress = zeros(2,2,n);
    exactStress(1,1,:) = reshape(4*(lambda+mu)*x,1,1,n);
    exactStress(2,2,:) = reshape(4*(lambda+mu)*x,1,1,n);
    exactStress(1,2,:) = reshape(4*mu*y,1,1,n);
    exactStress(2,1,:) = reshape(4*mu*y,1,1,n);

    displacementError = U-Uexact;
    strainError = strain-exactStrain;
    stressError = stress-exactStress;
    displacementMagnitude = sqrt(sum(displacementError.^2,2));
    strainMagnitude = tensor_nodal_norm(strainError);
    stressMagnitude = tensor_nodal_norm(stressError);
    sol.errors.displacement = error_statistics(displacementMagnitude, ...
        sqrt(sum(Uexact.^2,2)),mesh.vertexWeights,F,mesh.isBoundary);
    sol.errors.strain = error_statistics(strainMagnitude, ...
        tensor_nodal_norm(exactStrain),mesh.vertexWeights,F,mesh.isBoundary);
    sol.errors.stress = error_statistics(stressMagnitude, ...
        tensor_nodal_norm(exactStress),mesh.vertexWeights,F,mesh.isBoundary);

    residual = K*u-b;
    uExactVector = reshape(Uexact.',[],1);
    exactResidual = K*uExactVector-b;
    sol.residuals.solver = residual_statistics(residual,dofF,F,rhs,n);
    sol.residuals.analytical = residual_statistics(exactResidual,dofF,F,rhs,n);
    GfromExact = reconstruct_nodal_gradients(P,Uexact,mesh.adjacency);
    sol.exactSampleReconstructionGradient = GfromExact;
    sol.exactSampleReconstructionError = tensor_nodal_norm(GfromExact-Gexact);

    gradientError = G-Gexact;
    e11 = reshape(gradientError(1,1,:),n,1);
    e12 = reshape(gradientError(1,2,:),n,1);
    e21 = reshape(gradientError(2,1,:),n,1);
    e22 = reshape(gradientError(2,2,:),n,1);
    diagonalError = sqrt(e11.^2+e22.^2);
    offDiagonalError = sqrt(e12.^2+e21.^2);
    fullGradientError = sqrt(diagonalError.^2+offDiagonalError.^2);

    % Compare the diagonal and off-diagonal parts with their exact values.
    exactDiagonalMagnitude = sqrt(reshape(Gexact(1,1,:),n,1).^2 + ...
                                  reshape(Gexact(2,2,:),n,1).^2);
    exactOffDiagonalMagnitude = sqrt(reshape(Gexact(1,2,:),n,1).^2 + ...
                                     reshape(Gexact(2,1,:),n,1).^2);
    sol.errors.gradient = error_statistics(fullGradientError, ...
        tensor_nodal_norm(Gexact),mesh.vertexWeights,F,mesh.isBoundary);
    sol.errors.gradientDiagonal = error_statistics(diagonalError, ...
        exactDiagonalMagnitude,mesh.vertexWeights,F,mesh.isBoundary);
    sol.errors.gradientOffDiagonal = error_statistics(offDiagonalError, ...
        exactOffDiagonalMagnitude,mesh.vertexWeights,F,mesh.isBoundary);

    % Keep both shear derivatives separate: G_h12 and G_h21 need not agree.
    componentErrors = [e11,e22,e12,e21];
    exactComponents = [reshape(Gexact(1,1,:),n,1), ...
                       reshape(Gexact(2,2,:),n,1), ...
                       reshape(Gexact(1,2,:),n,1), ...
                       reshape(Gexact(2,1,:),n,1)];
    componentMean = mean(abs(componentErrors),1).';
    componentMax = max(abs(componentErrors),[],1).';
    componentL2 = sqrt(sum(bsxfun(@times,mesh.vertexWeights, ...
                                componentErrors.^2),1)).';
    exactComponentL2 = sqrt(sum(bsxfun(@times,mesh.vertexWeights, ...
                                     exactComponents.^2),1)).';
    componentRelativeL2 = zeros(4,1);
    for component = 1:4
        componentRelativeL2(component) = safe_ratio( ...
            componentL2(component),exactComponentL2(component));
    end
    componentLabels = {'G11 = du1/dx';'G22 = du2/dy'; ...
                       'G12 = du1/dy';'G21 = du2/dx'};
    sol.gradientComponentOrder = {'G11','G22','G12','G21'};
    sol.gradientComponentErrors = componentErrors;
    sol.gradientComponentErrorTable = table(componentLabels,componentMean, ...
        componentMax,componentL2,componentRelativeL2, ...
        'VariableNames',{'Component','MeanAbsolute','MaxAbsolute', ...
                         'WeightedL2','RelativeL2'});

    groupNames = {'Diagonal';'Off-diagonal';'Full gradient'};
    groupFields = {'gradientDiagonal','gradientOffDiagonal','gradient'};
    groupValues = zeros(3,4);
    for group = 1:3
        stats = sol.errors.(groupFields{group});
        groupValues(group,:) = [stats.mean,stats.maximum, ...
            stats.weightedL2,stats.relativeWeightedL2];
    end
    sol.gradientGroupedErrorTable = table(groupNames,groupValues(:,1), ...
        groupValues(:,2),groupValues(:,3),groupValues(:,4), ...
        'VariableNames',{'GradientPart','MeanNodal','MaxNodal', ...
                         'WeightedL2','RelativeL2'});

    sol.K = K;
    sol.load = b;
    sol.bodyForce = f;
    sol.boundaryTractionLoad = boundaryLoad;
    sol.analyticalCircleTraction = ...
        [4*(lambda+mu)*x.^2+4*mu*y.^2,4*(lambda+2*mu)*x.*y];
    sol.dirichletData = Uexact;
    sol.nodalDisplacement = U;
    sol.exactDisplacement = Uexact;
    sol.analyticalDeformation = P+Uexact;
    sol.discreteDeformation = P+U;
    sol.discreteStrain = strain;
    sol.exactStrain = exactStrain;
    sol.strainError = strainError;
    sol.nodalStrainError = strainMagnitude;
    sol.discreteStress = stress;
    sol.exactStress = exactStress;
    sol.stressError = stressError;
    sol.nodalStressError = stressMagnitude;
    sol.nodalError = displacementError;
    sol.nodalErrorMagnitude = displacementMagnitude;
    sol.nodalErrorMagnitudeMean = sol.errors.displacement.mean;
    sol.maxNodalError = sol.errors.displacement.maximum;
    sol.maxNodalStrainError = sol.errors.strain.maximum;
    sol.maxNodalStressError = sol.errors.stress.maximum;
    sol.weightedL2Error = sol.errors.displacement.weightedL2;
    sol.relativeWeightedL2Error = sol.errors.displacement.relativeWeightedL2;
    sol.dirichletError = max(displacementMagnitude(D));
    sol.freeResidual = sol.residuals.solver.relativeNorm2;
    sol.dirichletVertices = D;
    sol.freeVertices = F;
    sol.dirichletConstraintError = max(displacementMagnitude(D));
    sol.prescribedDirichletDeformation = P(D,:)+Uexact(D,:);
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
        [mean(abs(e11)),mean(abs(e12));mean(abs(e21)),mean(abs(e22))];
end


function norms = tensor_nodal_norm(tensor)
% Full Frobenius norm, not the three-entry norm of a Voigt vector.
    norms = sqrt(sum(reshape(tensor,4,[]).^2,1)).';
end


function stats = error_statistics(errors,exactMagnitudes,V,F,isBoundary)
    stats.nodal = errors;
    stats.maximum = max(errors);
    stats.maximumFree = max(errors(F));
    stats.maximumInterior = max(errors(~isBoundary));
    stats.maximumBoundary = max(errors(isBoundary));
    stats.mean = mean(errors);
    stats.weightedL2 = sqrt(sum(V.*errors.^2));
    stats.exactWeightedL2 = sqrt(sum(V.*exactMagnitudes.^2));
    stats.relativeWeightedL2 = safe_ratio(stats.weightedL2,stats.exactWeightedL2);
end


function stats = residual_statistics(residual,dofF,F,rhs,n)
    nodal = reshape(residual,2,n).';
    free = residual(dofF);
    stats.freeVector = free;
    stats.freeVertices = F;
    stats.freeNodalVectors = nodal(F,:);
    stats.freeNodalMagnitudes = sqrt(sum(nodal(F,:).^2,2));
    stats.nodalMagnitudes = NaN(n,1); % Dirichlet reactions are not free residuals.
    stats.nodalMagnitudes(F) = stats.freeNodalMagnitudes;
    stats.norm2 = norm(free);
    stats.relativeNorm2 = safe_ratio(stats.norm2,norm(rhs));
    stats.maxNodal = max(stats.freeNodalMagnitudes);
end


function ratio = safe_ratio(numerator,denominator)
    if denominator > 0
        ratio = numerator/denominator;
    elseif numerator == 0
        ratio = 0;
    else
        ratio = NaN;
    end
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
    fprintf('%s: N actual = %d, triangles = %d, h = %.6e\n', ...
        name,size(mesh.points,1),size(mesh.triangles,1),mesh.h);
    if isfield(mesh,'regularGridSpacing')
        fprintf('  Cartesian spacing = %.6e (boundary transition included).\n', ...
            mesh.regularGridSpacing);
    end
    fprintf('  Dirichlet vertices = %d, free vertices = %d\n', ...
        numel(sol.dirichletVertices),numel(sol.freeVertices));
    fprintf('  Field          weighted L2      relative L2      max nodal        mean nodal\n');
    fields = {'displacement','strain','stress'};
    for k = 1:numel(fields)
        e = sol.errors.(fields{k});
        fprintf('  %-13s  %.6e     %.6e     %.6e     %.6e\n', ...
            fields{k},e.weightedL2,e.relativeWeightedL2,e.maximum,e.mean);
    end
    r = sol.residuals.solver;
    a = sol.residuals.analytical;
    fprintf('  Solver free residual:     ||r_F||2=%.6e, relative=%.6e, max node=%.6e\n', ...
        r.norm2,r.relativeNorm2,r.maxNodal);
    fprintf('  Analytical free residual: ||r_F||2=%.6e, relative=%.6e, max node=%.6e\n', ...
        a.norm2,a.relativeNorm2,a.maxNodal);
    fprintf('  Dirichlet constraint error = %.6e\n',sol.dirichletConstraintError);
    fprintf('  Max gradient reconstruction error on exact samples = %.6e\n', ...
        max(sol.exactSampleReconstructionError));
    fprintf('  Gradient comparison: G_h - grad(u_exact), G_ij = d_j u_i.\n');
    fprintf('  Diagonal components: G11, G22; off-diagonal components: G12, G21.\n');
    disp(sol.gradientComponentErrorTable);
    fprintf('  Grouped nodal gradient-error norms:\n');
    disp(sol.gradientGroupedErrorTable);
end


function convergence = build_convergence(levels,family)
    m = numel(levels);
    targets = zeros(m,1);
    N = zeros(m,1);
    h = zeros(m,1);
    errors = zeros(m,6);
    residuals = zeros(m,6);
    gradientGroups = zeros(m,6);
    for k = 1:m
        mesh = levels(k).(family).mesh;
        sol = levels(k).(family).solution;
        targets(k) = levels(k).targetNodeCount;
        N(k) = size(mesh.points,1);
        h(k) = mesh.h;
        errors(k,:) = [sol.errors.displacement.weightedL2, ...
            sol.errors.displacement.maximum,sol.errors.strain.weightedL2, ...
            sol.errors.strain.maximum,sol.errors.stress.weightedL2, ...
            sol.errors.stress.maximum];
        r = sol.residuals.solver;
        a = sol.residuals.analytical;
        residuals(k,:) = [r.norm2,r.relativeNorm2,r.maxNodal, ...
                         a.norm2,a.relativeNorm2,a.maxNodal];
        diagonal = sol.errors.gradientDiagonal;
        offDiagonal = sol.errors.gradientOffDiagonal;
        gradientGroups(k,:) = [diagonal.mean,diagonal.maximum,diagonal.weightedL2, ...
            offDiagonal.mean,offDiagonal.maximum,offDiagonal.weightedL2];
    end
    orders = NaN(m,6);
    for k = 2:m
        if h(k) < h(k-1) && h(k) > 0
            for j = 1:6
                if errors(k-1,j) > 0 && errors(k,j) > 0 && ...
                        isfinite(errors(k-1,j)) && isfinite(errors(k,j))
                    orders(k,j) = log(errors(k-1,j)/errors(k,j))/ ...
                                  log(h(k-1)/h(k));
                end
            end
        else
            warning('DiskElasticity:NondecreasingH', ...
                '%s, levels %d/%d: h did not decrease; p_obs is NaN.',family,k-1,k);
        end
    end
    convergence.Ntarget = targets;
    convergence.Nactual = N;
    convergence.h = h;
    convergence.errors = errors;
    convergence.observedOrders = orders;
    convergence.errorColumns = {'u_L2','u_max','strain_L2','strain_max', ...
                                'stress_L2','stress_max'};
    convergence.errorTable = array2table([targets,N,h,errors], ...
        'VariableNames',[{'Ntarget','Nactual','h'},convergence.errorColumns]);
    convergence.orderTable = array2table([N,h,orders], ...
        'VariableNames',{'Nactual','h','p_u_L2','p_u_max','p_strain_L2', ...
                         'p_strain_max','p_stress_L2','p_stress_max'});
    convergence.residualTable = array2table([N,h,residuals], ...
        'VariableNames',{'Nactual','h','solver_norm2','solver_relative', ...
                         'solver_max_node','exact_norm2','exact_relative', ...
                         'exact_max_node'});
    convergence.gradientErrorTable = array2table([targets,N,h,gradientGroups], ...
        'VariableNames',{'Ntarget','Nactual','h','diag_mean','diag_max','diag_L2', ...
                         'offdiag_mean','offdiag_max','offdiag_L2'});
end


function print_convergence(name,convergence)
    fprintf('\n%s: weighted L2 and maximum nodal errors\n',name);
    disp(convergence.errorTable);
    fprintf('%s: p_obs between consecutive rows (first row NaN)\n',name);
    disp(convergence.orderTable);
    fprintf('%s: equilibrium residuals restricted to free vertices\n',name);
    disp(convergence.residualTable);
    fprintf('%s: diagonal and off-diagonal gradient errors at every level\n',name);
    disp(convergence.gradientErrorTable);
end


function plot_field_comparisons(meshRand,solRand,meshRegular,solRegular)
% Each figure compares exact norm, computed norm, and norm of the full error.
% The error panel is NOT the difference between the two displayed norms.
    meshes = {meshRand,meshRegular};
    solutions = {solRand,solRegular};
    names = {'Random','Regular'};
    labels = {'Displacement','Strain','Stress'};
    for field = 1:3
        exact = cell(2,1);
        computed = cell(2,1);
        errors = cell(2,1);
        for row = 1:2
            sol = solutions{row};
            if field == 1
                exact{row} = sqrt(sum(sol.exactDisplacement.^2,2));
                computed{row} = sqrt(sum(sol.nodalDisplacement.^2,2));
                errors{row} = sol.nodalErrorMagnitude;
            elseif field == 2
                exact{row} = tensor_nodal_norm(sol.exactStrain);
                computed{row} = tensor_nodal_norm(sol.discreteStrain);
                errors{row} = sol.nodalStrainError;
            else
                exact{row} = tensor_nodal_norm(sol.exactStress);
                computed{row} = tensor_nodal_norm(sol.discreteStress);
                errors{row} = sol.nodalStressError;
            end
        end
        fieldLimit = max([exact{1};exact{2};computed{1};computed{2}]);
        errorLimit = max([errors{1};errors{2}]);
        create_large_figure([labels{field} ': analytical and discrete comparison']);
        for row = 1:2
            subplot(2,3,3*(row-1)+1);
            draw_gradient_error(meshes{row},exact{row},max(fieldLimit,eps));
            title([names{row} ': analytical ' lower(labels{field}) ' norm']);
            subplot(2,3,3*(row-1)+2);
            draw_gradient_error(meshes{row},computed{row},max(fieldLimit,eps));
            title([names{row} ': discrete ' lower(labels{field}) ' norm']);
            subplot(2,3,3*(row-1)+3);
            draw_gradient_error(meshes{row},errors{row},max(errorLimit,eps));
            title({[names{row} ': full error norm'], ...
                sprintf('Maximum = %.3e',max(errors{row}))});
        end
    end
end


function plot_free_residuals(meshRand,solRand,meshRegular,solRegular)
    meshes = {meshRand,meshRegular};
    solutions = {solRand,solRegular};
    names = {'Random','Regular'};
    types = {'solver','analytical'};
    create_large_figure('Free-vertex equilibrium residuals (nodal force norms)');
    for column = 1:2
        r1 = solutions{1}.residuals.(types{column});
        r2 = solutions{2}.residuals.(types{column});
        limit = max([r1.maxNodal;r2.maxNodal;eps]);
        for row = 1:2
            mesh = meshes{row};
            r = solutions{row}.residuals.(types{column});
            F = r.freeVertices;
            subplot(2,2,2*(row-1)+column);
            scatter(mesh.points(F,1),mesh.points(F,2),24, ...
                    r.freeNodalMagnitudes,'filled');
            axis equal tight; colorbar; caxis([0,limit]);
            xlabel('x'); ylabel('y');
            title({[names{row} ': ' types{column} ' free residual'], ...
                sprintf('Max = %.2e; relative = %.2e',r.maxNodal,r.relativeNorm2)});
        end
    end
end


function plot_convergence(convergence)
    create_large_figure('Observed convergence: actual mesh size h');
    labels = {'Displacement L2','Displacement maximum', ...
              'Strain L2','Strain maximum','Stress L2','Stress maximum'};
    for k = 1:6
        subplot(3,2,k);
        loglog(convergence.random.h,convergence.random.errors(:,k), ...
               'o-','LineWidth',1.3); hold on;
        loglog(convergence.regular.h,convergence.regular.errors(:,k), ...
               's-','LineWidth',1.3); hold off;
        grid on; xlabel('h = maximum edge length'); ylabel('Error');
        title(labels{k});
        legend('Random','Regular','Location','best');
    end
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
    title('Regular point distribution');
    xlabel('x'); ylabel('y');
    legend('Delaunay edges','\Gamma_D','\Gamma_N','Location','southoutside');
end


function plot_deformations(meshRand,solRand,meshEq,solEq)
%PLOT_DEFORMATIONS Group six point-cloud views in one large 2-by-3 figure.
% Rows are random/regular. Columns are discrete deformation,
% original disk, and analytical deformation alone. All views retain the
% same reference-node colors and use Phi = X + U with displacement scale 1.
% The color key uses the original angle theta_0 and radius r_0, never the
% deformed coordinates. Six configuration panels remain in a 2-by-3 layout.
% Each panel fits its own points, with no shared limits or added padding.
% Separate header bands keep captions clear of axes and tick labels.

    meshes = {meshRand,meshEq};
    solutions = {solRand,solEq};
    names = {'Random','Regular'};

    create_large_figure('Six views: original disks and deformations');
    deformationFigure = gcf;
    columnLeft = [0.035,0.355,0.675];
    panelWidth = 0.29;
    panelHeight = 0.335;
    for k = 1:2
        mesh = meshes{k};
        sol = solutions{k};
        colors = reference_node_colors(mesh.points);
        % Leave a bottom strip for the two-dimensional color key.
        panelBottom = 0.57 - 0.40*(k-1);
        configurations = {sol.discreteDeformation,mesh.points, ...
                          sol.analyticalDeformation};
        captions = {'Discrete deformation','Original disk', ...
                    'Analytical deformation'};
        for column = 1:3
            outerPosition = [columnLeft(column),panelBottom, ...
                             panelWidth,panelHeight];
            panelAxes = axes('Parent',deformationFigure,'Units','normalized');
            draw_colored_point_cloud(configurations{column},colors,mesh.isBoundary);
            % OuterPosition includes ticks and their labels. Use the R2017b
            % property so these decorations stay inside the allocated slot.
            set(panelAxes,'ActivePositionProperty','outerposition', ...
                'OuterPosition',outerPosition);
            headerPosition = [columnLeft(column), ...
                panelBottom+panelHeight+0.010,panelWidth,0.040];
            annotation(deformationFigure,'textbox',headerPosition, ...
                'String',[names{k} ': ' captions{column}], ...
                'Units','normalized','FitBoxToText','off','Margin',2, ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'EdgeColor','none','FontSize',11,'FontWeight','bold', ...
                'Interpreter','none');
        end
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


function draw_colored_point_cloud(points,colors,isBoundary)
%DRAW_COLORED_POINT_CLOUD Plot a single configuration with reference colors.
% Draw interiors first so that boundary nodes remain visible in both meshes.

    markerArea = max(8,min(32,32000/size(points,1)));
    interior = ~isBoundary;
    scatter(points(interior,1),points(interior,2),markerArea, ...
        colors(interior,:),'o','filled', ...
        'MarkerEdgeColor','none','Clipping','off');
    hold on;
    scatter(points(isBoundary,1),points(isBoundary,2),1.5*markerArea, ...
        colors(isBoundary,:),'o','filled', ...
        'MarkerEdgeColor',[0.20 0.20 0.20],'LineWidth',0.5,'Clipping','off');
    hold off;
    % Fit this configuration alone while preserving geometric proportions.
    % Clipping is off so markers at the extrema remain completely visible.
    axis equal tight;
    set(gca,'FontSize',10);
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
        title({['Regular: ' labels{r}], ...
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
