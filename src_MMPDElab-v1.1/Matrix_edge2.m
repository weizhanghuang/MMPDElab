function [E,Einv,detE,xK,h] = Matrix_edge2(X,tri,j1,k1)
%
% usage: [E,Einv,detE,xK,h] = Matrix_edge2(X,tri,j1,k1)
%
% this function is similar to Matrix_edge(X,tri), i.e., computes
% the edge matrix, its inverse and determinant, and barycenters of the elements
% for mesh (X, tri), the k1-th component of the coordinates of j1-th vertex
% of all elements perturbed. this function is needed in the computation
% of Jacobian matrix related to mesh velocities.
%
% X:    coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  connectivity of the mesh, of size N-by-(d+1).
% E:    (output) edge matrix, E = [(X_2-X_1)_1, ..., (X_2-X_1)_d, ...,
%       (X_{d+1}-X_1)_1, ..., (X_{d+1}-X_1)_d], of size N-by-d*d.
% Einv: (output) inverse of E, of size N-by-d*d.
% detE: (output) determinant of the edge matrix, of size N-by-1.
% xK:   (output) barycenters of the elements, of size N-by-d.
% h:    (output) pertubation in the k1-th component of the coordinates of the
%       j1-th vertex of all elements, of size N-by-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

[N,d] = size(tri);
d = d - 1;

EE = zeros(N,d,d+1);
for j=1:(d+1)
   EE(:,:,j) = X(tri(:,j),:);
end

h = (EE(:,k1,j1)+max(abs(EE(:,k1,j1)),1)*sqrt(eps))-EE(:,k1,j1);
ss = (sign(EE(:,k1,j1))>= 0);
h = (ss-(~ss)).*abs(h);

EE(:,k1,j1) = EE(:,k1,j1) + h; 
 
xK = zeros(N,d);
for j=1:(d+1)
   xK = xK + EE(:,:,j);
end
xK = xK/(d+1);
          
for j=2:(d+1)
   EE(:,:,j) = EE(:,:,j) - EE(:,:,1);
end
E = zeros(N,d*d);
for j=2:(d+1)
for k=1:d
   E(:,d*(j-2)+k) = EE(:,k,j);
end
end
detE = Matrix_det(E);
Einv = Matrix_inv(E,detE);

% end of Matrix_edge2()
