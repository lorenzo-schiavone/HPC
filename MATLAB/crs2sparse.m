% function A = crs2sparse(nnz, iat, ja, coef)
% % essentially i have to convert iat in rows and add one to each ja
% rows = zeros(nnz, 1);
% nrows = length(iat)-1;
% ncols = nrows;
% for i=1:nrows
%     for j=iat(i):(iat(i+1)-1)
%     rows(j+1) = i;
%     end
% end
% A = sparse(rows,ja+1, coef, nrows, ncols, nnz);
% end 

function A = crs2sparse(nnz, iat, ja, coef)
    % Convert CRS format (0-based indexing) to MATLAB sparse matrix
    rows = zeros(nnz, 1);
    nrows = length(iat) - 1;
    ncols = nrows;

    for i = 1:nrows
        for j = iat(i)+1 : iat(i+1) % +1 for MATLAB 1-based indexing
            rows(j) = i;
        end
    end

    A = sparse(rows, ja + 1, coef, nrows, ncols, nnz); % ja+1 to move to 1-based
end
