clc; close all; clear;
elements = load("Cubo_591.tetra");
nodes = load("Cubo_591.coor");


D = [1,1,1];
v = [1,1,1];

tot_el = size(elements,1);
tot_nodes = size(nodes, 1);
volumes = zeros(tot_el, 1);
tic
%% computing areas of each element
parfor i = 1:tot_el
    matrix = zeros(4,4);
    for j =1:4
        matrix(j,:)= [1 nodes(elements(i, j), :)];
    end
    volumes(i,1)=det(matrix)/6;
end


I = cell(tot_el,1); J = cell(tot_el,1); 
Hval = cell(tot_el,1);
Pval = cell(tot_el,1); Bval = cell(tot_el,1);

for kk = 1
    el_nodes = elements(kk, :);
    vol = volumes(kk);

    a = zeros(4,1);
    b = zeros(4,1);
    c = zeros(4,1);
    d = zeros(4,1);
    for i=1:4
        j = mod(i,4)+1; k = mod(i+1,4)+1; m = mod(i+2,4)+1;
        fprintf("%u %u %u\n", el_nodes(j),el_nodes(k),el_nodes(m));
        nj = nodes(el_nodes(j),:);
        nk = nodes(el_nodes(k),:);
        nm = nodes(el_nodes(m),:);
        M = [nj;nk;nm];
        M 
        a(i) = det(M);
        b(i) = -det([1, nj(2), nj(3);
             1, nk(2), nk(3);
             1, nm(2), nm(3)]);
        c(i) = det([1, nj([1,3]); ...
                     1, nk([1,3]); ...
                     1, nm([1,3]);
                     ]);
        d(i) = - det([1, nj(1:2); ...
                     1, nk(1:2); ...
                     1, nm(1:2);
                     ]);
    end
    vol
    a
    b
    c
    d

    ...
    Hloc = (D(1)*(b*b') + D(2)*(c*c') + D(3)*(d*d'))/(36*abs(vol));
    Ploc = abs(vol)/20 * (ones(4) + eye(4));
    bloc = (v(1)*b+v(2)*c+v(3)*d)/(24 *sign(vol));
    Bloc = [bloc, bloc, bloc, bloc];
    [ii, jj] = meshgrid(el_nodes, el_nodes);
    I{kk} = ii(:); J{kk} = jj(:);
    Hval{kk} = Hloc(:);
    Pval{kk} = Ploc(:);
    Bval{kk} = Bloc(:);
end

H = sparse(cell2mat(I), cell2mat(J), cell2mat(Hval), tot_nodes, tot_nodes);
P = sparse(cell2mat(I), cell2mat(J), cell2mat(Pval), tot_nodes, tot_nodes);
B = sparse(cell2mat(I), cell2mat(J), cell2mat(Bval), tot_nodes, tot_nodes);


load("coefB.txt")
load("coefH.txt")
load("coefP.txt")
load("iat.txt")
load("ja.txt")
load("nnz.txt")

B1 = crs2sparse(nnz,iat,ja, coefB);
P1 = crs2sparse(nnz,iat,ja, coefP);
H1 = crs2sparse(nnz,iat,ja, coefH);

disp(sum(sum(abs(H1-H)))/nnz)
disp(sum(sum(abs(P1-P)))/nnz)
disp(sum(sum(abs(B1-B)))/nnz)