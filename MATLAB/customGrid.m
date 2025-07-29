function [nodes, elements] = generateCubeTetraMesh(L, nx, ny, nz)
% generateCubeTetraMesh: create a structured tetrahedral mesh of a cube.
%
% INPUTS:
%   L      - Length of the cube edge (cube = [0,L]^3)
%   nx,ny,nz - Number of subdivisions in x, y, z directions
%
% OUTPUTS:
%   nodes     - N×3 matrix of node coordinates
%   elements  - M×4 matrix of tetrahedron vertex indices

% Create grid of points
[x, y, z] = ndgrid(linspace(0, L, nx+1), ...
                   linspace(0, L, ny+1), ...
                   linspace(0, L, nz+1));

nodes = [x(:), y(:), z(:)];
nodeID = reshape(1:size(nodes,1), size(x));  % map 3D indices to global IDs

elements = [];

% Loop through each cube cell and divide into tetrahedra
for i = 1:nx
    for j = 1:ny
        for k = 1:nz
            % Indices of the 8 corners of the current cube
            n000 = nodeID(i,   j,   k);
            n100 = nodeID(i+1, j,   k);
            n010 = nodeID(i,   j+1, k);
            n110 = nodeID(i+1, j+1, k);
            n001 = nodeID(i,   j,   k+1);
            n101 = nodeID(i+1, j,   k+1);
            n011 = nodeID(i,   j+1, k+1);
            n111 = nodeID(i+1, j+1, k+1);

            % Divide cube into 6 tetrahedra (standard pattern)
            elements = [elements;
                n000, n001, n011, n111;
                n000, n011, n010, n111;
                n000, n010, n110, n111;
                n000, n110, n100, n111;
                n000, n100, n101, n111;
                n000, n101, n001, n111];
        end
    end
end
end


L = 1;        % unit cube
nx = 20; ny = 20; nz = 20;  % number of divisions per direction

[nodes, elements] = generateCubeTetraMesh(L, nx, ny, nz);

% % Visualize
% tetramesh(elements, nodes);
% axis equal; view(3);
% title('Structured Tetrahedral Mesh of a Cube');


function saveMeshCoorTetra(nodes, elements, basename)
% Save FEM mesh in .coor and .tetra format
%
% Inputs:
%   nodes    - N×3 matrix of node coordinates
%   elements - M×4 matrix of tetrahedral element node indices
%   basename - base filename (e.g., 'CubeMesh')
%
% Output:
%   Writes basename.coor and basename.tetra

% Number of nodes and elements
num_nodes = size(nodes, 1);
num_elements = size(elements, 1);

% Save .coor file
coor_file = fopen([basename, '.coor'], 'w');
for i = 1:num_nodes
    fprintf(coor_file, '%.15g %.15g %.15g\n', nodes(i,1), nodes(i,2), nodes(i,3));
end
fclose(coor_file);

% Save .tetra file (region marker set to 1)
tetra_file = fopen([basename, '.tetra'], 'w');
for i = 1:num_elements
    fprintf(tetra_file, '%d %d %d %d 1\n', elements(i,1), elements(i,2), elements(i,3), elements(i,4));
end
fclose(tetra_file);

fprintf('Mesh saved to %s.coor and %s.tetra\n', basename, basename);
end


saveMeshCoorTetra(nodes, elements, sprintf('CubeMesh_%d', nx));