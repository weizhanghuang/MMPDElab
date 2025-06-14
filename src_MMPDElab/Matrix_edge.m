function [E,Einv,detE,xK] = Matrix_edge(X,tri)
%
% usage: [E,Einv,detE,xK] = Matrix_edge(X,tri)
%
% this function computes the edge matrix, its inverse and determinant, and
% barycenters of the elements for mesh (X, tri).
%
% X:    coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  connectivity of the mesh, of size N-by-(d+1).
% E:    (output) edge matrix, E = [(X_2-X_1)_1, ..., (X_2-X_1)_d, ...,
%       (X_{d+1}-X_1)_1, ..., (X_{d+1}-X_1)_d], of size N-by-d*d.
% Einv: (output) inverse of E, of size N-by-d*d.
% detE: (output) determinant of the edge matrix, of size N-by-1.
% xK:   (output) barycenters of the elements, of size N-by-d.
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
Einv = Matrix_inv(E,detE);

xK = zeros(N,d);
for j=1:(d+1)
   xK = xK + X(tri(:,j),:);
end
xK = xK/(d+1);

% end of Matrix_edge()
