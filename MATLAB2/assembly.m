clc; close all; clear;
elements = load("Cubo_4820_rcm.tetra");
nodes = load("Cubo_4820_rcm.coor");

% % to have positive volumes
% el_temp = elements(:,:);
% elements(:,1) = el_temp(:,2);
% elements(:,2) = el_temp(:,1);

% elements = elements(:,2:end-1);
% nodes = nodes(:, 2:end);

D = 1e-1 * [4,1,1];
v = [1,1,2];
dt = 1;
nsteps = 1;
tol = 1e-9;
maxit = 1;
repeat = 2000;


tot_el = size(elements,1);
tot_nodes = size(nodes, 1);
volumes = zeros(tot_el, 1);

tic
% computing volumes of each element
parfor i = 1:tot_el
    matrix = [ 1 nodes(elements(i, 1), :);
               1 nodes(elements(i, 2), :);
               1 nodes(elements(i, 3), :);
               1 nodes(elements(i, 4), :);
              ];
    volumes(i,1)=det(matrix)/6;
end

H = sparse(tot_nodes, tot_nodes);
P = sparse(tot_nodes, tot_nodes);
B = sparse(tot_nodes, tot_nodes);

tic 
for kk = 1:tot_el
    el_nodes = elements(kk, 1:4);
    vol = volumes(kk);
    a = zeros(4,1);
    b = zeros(4,1);
    c = zeros(4,1);
    d = zeros(4,1);
    for i=1:4
        %j = mod(i,4)+1; k = mod(i+1,4)+1; m = mod(i+2,4)+1; WRONG
        idxs = setdiff(1:4,i); j = idxs(1); k = idxs(2); m = idxs(3);
        sig = (-1)^(i+1);
        nj = nodes(el_nodes(j),:);
        nk = nodes(el_nodes(k),:);
        nm = nodes(el_nodes(m),:);

        a(i) = sig * det([nj;nk;nm]);
        b(i) = sig * (-1) * det([1, nj(2), nj(3);
                     1, nk(2), nk(3);
                     1, nm(2), nm(3); 
                     ]);
        c(i) = sig * det([1, nj([1,3]); 
                    1, nk([1,3]);
                    1, nm([1,3]);
                    ]);
        d(i) = sig * (-1) * det([1, nj(1:2);
                      1, nk(1:2);
                      1, nm(1:2);
                      ]);
    end

    % if kk==77
    %     disp("a: ");
    %     disp(a);
    % 
    %     disp("b: ");
    %     disp(b);
    % 
    %     disp("c: ");
    %     disp(c);
    % 
    %     disp("d: ");
    %     disp(d);
    % end

    % Xcord = [ 1 nodes(elements(kk, 1), :);
    %            1 nodes(elements(kk, 2), :);
    %            1 nodes(elements(kk, 3), :);
    %            1 nodes(elements(kk, 4), :);
    %           ];
    % 
    % inv_Xcord = 6 * vol * inv(Xcord);
    % a_new = inv_Xcord(1,:)';
    % b_new = inv_Xcord(2,:)';
    % c_new = inv_Xcord(3,:)';
    % d_new = inv_Xcord(4,:)';
    
    Hloc = ( D(1)*(b*b') + D(2)*(c*c') + D(3)*(d*d')) / (36*abs(vol));
    Ploc = abs(vol)/20 * (ones(4) + eye(4));
    Bloc = sign(vol)/24 * repmat(v * [b,c,d]', 4,1); % sign(vol)
    
    for i = 1:4
        row = el_nodes(i);
        for j=1:4
            col = el_nodes(j);
            H(row,col) = H(row,col) + Hloc(i,j);
            P(row,col) = P(row,col) + Ploc(i,j);
            B(row,col) = B(row,col) + Bloc(i,j);
        end
    end
end
toc
disp("Assemebly completed!");

% % CHECK IF EQUAL TO CPP
% load("coefB.txt")
% load("coefH.txt")
% load("coefP.txt")
% load("iat.txt")
% load("ja.txt")
% load("nnz.txt")
% 
% B1 = crs2sparse(nnz,iat,ja, coefB);
% P1 = crs2sparse(nnz,iat,ja, coefP);
% H1 = crs2sparse(nnz,iat,ja, coefH);
% 
% disp(sum(sum(abs(H1-H)))/nnz)
% disp(sum(sum(abs(P1-P)))/nnz)
% disp(sum(sum(abs(B1-B)))/nnz)

A_orig = B + H + P/dt;  

figure(200)
spy(A_orig)

% IC or BC
% r=.1;
% bc_nodes = find((abs(nodes(:,1)-.5) < r + 1e-10) & (abs(nodes(:,2)-.5) < r + 1e-10) & (abs(nodes(:,3)-.5) < r + 1e-10));
bc_nodes = find( (nodes(:,1)<1e-10) & (nodes(:,2)<1e-10) & (nodes(:,3)<.3 + 1e-10));
tot_bc_nodes = length(bc_nodes);
bc_val = ones(tot_bc_nodes, 1);

% % Precompute boundary adjustment term (A_orig * u_bc)
boundary_adjust = A_orig(:, bc_nodes) * bc_val;

% % Create modified system matrix (with boundary conditions applied)
A_mod = A_orig;
A_mod(bc_nodes, :) = 0;
A_mod(:, bc_nodes) = 0;
A_mod(bc_nodes, bc_nodes) = speye(tot_bc_nodes);

% Preconditioner 
%[L, U] = ilu(A_mod, struct('type', 'ilutp', 'droptol', 1e-3, 'udiag',1));
% Jacobi preconditioner (as in C)
Dinv = 1 ./ diag(A_mod);
M1 = @(x) Dinv .* x;

% Time-stepping loop
u_prev = zeros(tot_nodes, 1);  % Initial condition
u_prev(bc_nodes) = bc_val;

for step = 1:nsteps
    rhs = P/dt * u_prev - boundary_adjust;  
    rhs(bc_nodes) = bc_val;                 
    
    % Solve linear system
    
    u = gmres(A_mod, rhs, repeat , tol, maxit, M1);
    u_prev = u;  % Update for next time step
    
    % VISUALIZATION
    figure(step);
    trisurf(elements(:,1:4), nodes(:,1), nodes(:,2), nodes(:,3), u, 'FaceAlpha', .5);
    colorbar; %clim([0 1]); 
    shading interp;
    view(3); axis equal;
    drawnow; 
    xlabel("x");
    ylabel("y");
    zlabel("z");
    pause(0.1); 
end
disp("first 20 entries of the solution:")
disp(u(1:20))




