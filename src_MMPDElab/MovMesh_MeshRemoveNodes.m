function [X,tri] = MovMesh_MeshRemoveNodes(X1,tri1,ID)
%
% usage: [X,tri] = MovMesh_MeshRemoveNodes(X1,tri1,ID)
%
% this function removes the nodes (ID) from the existing mesh (X1,tri1).
%
% X1:   the coordinates of vertices of mesh 1, of size Nv1-by-d.
% tri1: the connectivity of mesh 1, of size N1-by-(d+1).
% ID:   IDs of the nodes to be removed.
% X:    (output) the coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  (output) the connectivity of the mesh, of size N-by-(d+1).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   ID1 = setdiff((1:size(X1,1))',ID);
   X = X1(ID1,:);
   d = size(X,2);
   ID2 = find(ismember(tri1(:,1),ID));
   for i=2:d+1
      ID2 = [ID2;find(ismember(tri1(:,i),ID))];
   end
   tri = tri1(setdiff((1:size(tri1,1))',ID2),:);
   II = zeros(size(X1,1),1);
   II(ID1) = (1:size(ID1,1))';
   tri = II(tri);
   
% end of MovMesh_MeshRemoveNodes()
