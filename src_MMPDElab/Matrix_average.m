function mK = Matrix_average(M,tri)
%
% usage: mK = Matrix_average(M,tri)
%
% this function computes the element-wise average of M defined
% on a triangular mesh with connectivity given by tri().
%
% M:    a vector of size Nv-by-m, where Nv is the number of vertices.
% tri:  connectivity of the mesh, of size N-by-(d+1).
% mK:   (output) a vector of size N-by-m, where N is the number of elements.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

[N,d] = size(tri);
d = d - 1;
m = size(M,2);

mK = zeros(N,m);
for j=1:(d+1)
   mK = mK + M(tri(:,j),:);
end
mK = mK/(d+1);

% end of Matrix_average()
