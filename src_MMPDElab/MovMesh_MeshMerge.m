function [X,tri] = MovMesh_MeshMerge(X1,tri1,X2,tri2)
%
% usage: [X,tri] = MovMesh_MeshMerge(X1,tri1,X2,tri2)
%
% this function merges two non-overlapping meshes (X1,tri1) and (X2,tri2).
% (they can share boundary nodes but non-overlapping.)
%
% X1:   the coordinates of vertices of mesh 1, of size Nv1-by-d.
% tri1: the connectivity of mesh 1, of size N1-by-(d+1).
% X2:   the coordinates of vertices of mesh 2, of size Nv2-by-d.
% tri2: the connectivity of mesh 2, of size N2-by-(d+1).
% X:    (output) the coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  (output) the connectivity of the mesh, of size N-by-(d+1).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   % merge
   X = [X1;X2];
   tri =[tri1;tri2+size(X1,1)];
   % removed repeated nodes
   [X,~,IC] = uniquetol(X,1e-8,'ByRows',true);
   tri = IC(tri);
   
% end of MovMesh_MeshMerge()
