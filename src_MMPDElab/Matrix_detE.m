function detE = Matrix_detE(X,tri)
%
% usage: detE = Matrix_detE(X,tri)
%
% this function computes the determinant of the edge matrices for mesh (X, tri).
%
% X:    coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  connectivity of the mesh, of size N-by-(d+1).
% detE: (output) determinant of the edge matrix, of size N-by-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

[N,d] = size(tri);
d = d - 1;

E = zeros(N,d*d);
for j=2:(d+1)
for k=1:d
   E(:,d*(j-2)+k) = X(tri(:,j),k)-X(tri(:,1),k);
end
end
detE = Matrix_det(E);

% end of Matrix_detE()
