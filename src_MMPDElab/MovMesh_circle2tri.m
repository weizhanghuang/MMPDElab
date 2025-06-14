function [X,tri] = MovMesh_circle2tri(jmax)
%
% usage: [X,tri] = MovMesh_circle2tri(jmax)
%
% this function generates a triangular mesh for the unit circle centered
% at (0,0) with radius 1.
%
% jmax: positive integer.
% X:    (output) coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  (output) the connectivity of the mesh, of size N-by-(d+1).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   [xx,yy] = ndgrid(linspace(-1,1,jmax),linspace(-1,1,jmax));
   X = [xx(:), yy(:)];
   X = X(X(:,1).^2+X(:,2).^2 < 1-1/jmax,:);
   jj = floor(3.2*jmax);
   X = [X; cos(linspace(0,2*pi*jj/(jj+1),jj)'), ...
           sin(linspace(0,2*pi*jj/(jj+1),jj)')];
   X = unique(X,'rows');
   TR = delaunayTriangulation(X);
   tri = TR.ConnectivityList;

% end of MovMesh_circle2tri()
