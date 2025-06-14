function bfE = MovMesh_bfAttachments(X,tri,tri_bf)
%
% usage: elm = MovMesh_bfAttachments(tri,tri_bf)
%
% this function identifies the elements attached to the given boundary facets.
%
% X:    the coordinates of vertices of the mesh of size Nv-by-d.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
% tri_bf: the boundary facets for all meshes, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% bfE:  (output) indices of the elements attached to the boundary facets,
%       of size Nbf-by-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2020  Weizhang Huang (whuang@ku.edu)

    [N,d] = size(tri);
    d = d-1;
    
    switch (d)
    case 1
       tri_bf1 = zeros(2*N,d);
       tri_bf1(1:N,:) = tri(:,1);
       tri_bf1(N+1:2*N,:) = tri(:,2);
       [Lia,Lib] = ismember(tri_bf1,tri_bf);
       [~,pp] = sort(Lib(Lia));
       ind = find(Lia);
       ind = ind(pp);
       bfE = mod(ind,N);
       bfE(bfE == 0) = N;
    case 2
%{
       tri_bf1 = zeros(3*N,d);
       tri_bf1(1:N,:) = [tri(:,1),tri(:,2)];
       tri_bf1(N+1:2*N,:) = [tri(:,1),tri(:,3)];
       tri_bf1(2*N+1:3*N,:) = [tri(:,2),tri(:,3)];
       tri_bf1 = sort(tri_bf1,2);
       [Lia,Lib] = ismember(tri_bf1,sort(tri_bf,2),'rows');
       [~,pp] = sort(Lib(Lia));
       ind = find(Lia);
       ind = ind(pp);
       bfE = mod(ind,N);
       bfE(bfE == 0) = N;
%}
       TR = triangulation(tri,X);
       bfElements = edgeAttachments(TR,tri_bf);
       bfE = cell2mat(bfElements);
       % bfE = cellfun(@(v) v(1),bfElements);
    case 3
       tri_bf1 = zeros(4*N,d);
       tri_bf1(1:N,:) = [tri(:,1),tri(:,2),tri(:,3)];
       tri_bf1(N+1:2*N,:) = [tri(:,1),tri(:,2),tri(:,4)];
       tri_bf1(2*N+1:3*N,:) = [tri(:,1),tri(:,3),tri(:,4)];
       tri_bf1(3*N+1:4*N,:) = [tri(:,2),tri(:,3),tri(:,4)];
       tri_bf1 = sort(tri_bf1,2);
       [Lia,Lib] = ismember(tri_bf1,sort(tri_bf,2),'rows');
       [~,pp] = sort(Lib(Lia));
       ind = find(Lia);
       ind = ind(pp);
       bfE = mod(ind,N);
       bfE(bfE == 0) = N;
%{
       TR = triangulation(tri,X);
       Elem1 = edgeAttachments(TR,tri_bf(:,1),tri_bf(:,2));
       bfElements = edgeAttachments(TR,tri_bf(:,2),tri_bf(:,3));
       bfElements = cellfun(@intersect,Elem1,bfElements,'UniformOutput',false);
       Elem1 = edgeAttachments(TR,tri_bf(:,1),tri_bf(:,3));
       bfElements = cellfun(@intersect,Elem1,bfElements,'UniformOutput',false);
       bfE = cell2mat(bfElements);
       % bfE = cellfun( @(v) v(1),bfElements);
%}
    end

% end of MovMesh_MovMesh_bfAttachments()
