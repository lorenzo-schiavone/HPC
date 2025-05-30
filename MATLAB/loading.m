load("coefB.txt")
load("coefH.txt")
load("coefP.txt")
load("iat.txt")
load("ja.txt")
load("nnz.txt")

B = crs2sparse(nnz,iat,ja, coefB);
P = crs2sparse(nnz,iat,ja, coefP);
H = crs2sparse(nnz,iat,ja, coefH);