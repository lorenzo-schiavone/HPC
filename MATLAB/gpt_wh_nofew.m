% Load mesh data
nodes = load('Cubo_591.coor');     % Nx3
elements = load('Cubo_591.tetra'); % Mx4

n_nodes = size(nodes, 1);
n_elem = size(elements, 1);

% Physical parameters
D = 1e-3;               % isotropic diffusion coefficient
v = [0, 0, 2];          % constant velocity field
dt = 0.1;              % time step
T_end = 1;              % final time
n_steps = T_end / dt;

% Allocate global matrices
H = sparse(n_nodes, n_nodes);  % stiffness (diffusion)
B = sparse(n_nodes, n_nodes);  % convection
P = sparse(n_nodes, n_nodes);  % capacity (mass)

% Assembly
for e = 1:n_elem
    % Node indices and coordinates
    idx = elements(e,:);
    Xe = nodes(idx, :);
    
    % Compute volume
    Ve = abs(det([ones(4,1), Xe])) / 6;
    
    % Compute gradients of basis functions
    A = [ones(4,1), Xe];
    coeffs = inv(A);  % columns: [a_i; b_i; c_i; d_i]
    b = coeffs(2, :)';
    c = coeffs(3, :)';
    d = coeffs(4, :)';
    
    % Local stiffness matrix (diffusion)
    He = (D/36) * (b*b' + c*c' + d*d') * abs(det(A));
    
    % Local mass matrix
    Pe = (Ve / 20) * (ones(4) + eye(4));
    
    % Local convection matrix
    vx = v(1); vy = v(2); vz = v(3);
    % Be = (1/24) * (vx * b' + vy * c' + vz * d') * ones(4,1)';
    % Be = Be .* Ve;
    Be = (1/24) * (vx * b + vy * c + vz * d) * ones(1,4);
    Be = Be * Ve;

    % Assemble into global matrices
    H(idx, idx) = H(idx, idx) + He;
    P(idx, idx) = P(idx, idx) + Pe;
    B(idx, idx) = B(idx, idx) + Be;
end

% Initial condition
c = zeros(n_nodes, 1);
% For example, a localized initial condition:
% c(vecnorm(nodes - center, 2, 2) < radius) = 1;
A_orig = B + H + P/dt;  % Store original A

% Dirichlet boundary (e.g., all boundary nodes to 0)
bc_nodes = find(((nodes(:,1) < 1e-10 & nodes(:,2) < 1+ 1e-10)| (nodes(:,1) < .3 + 1e-10 & nodes(:,2) < 1e-10) ) & nodes(:,3) > .99+ 1e-10);
tot_bc_nodes = length(bc_nodes);
bc_val = 10*ones(tot_bc_nodes, 1);
% bc_nodes = [1 2 3];
% tot_bc_nodes = 3;
% bc_val = [10;10;10];

% Precompute boundary adjustment term (A_orig * u_bc)
boundary_adjust = A_orig(:, bc_nodes) * bc_val;

% Create modified system matrix (with boundary conditions applied)
A_mod = A_orig;
A_mod(bc_nodes, :) = 0;
A_mod(:, bc_nodes) = 0;
A_mod(bc_nodes, bc_nodes) = speye(tot_bc_nodes);

% Time stepping (backward Euler)
A_sys = H + B + (1/dt)*P;
P_over_dt = (1/dt)*P;

for t = 1:n_steps
    rhs = P/dt * c - boundary_adjust;  % Includes time-stepping and BC terms
    rhs(bc_nodes) = bc_val;                 % Enforce BC values
    c_new = A_mod \ rhs;
    c = c_new;
    figure()
    trisurf(elements, nodes(:,1), nodes(:,2), nodes(:,3), c);
    axis equal; colorbar;
    title('Final concentration');
    view(3);
    pause(0.1)
end



