function V = MovMesh_freeBoundary_vertexNormal(X,tri,tri_bf)
%
% usage: V = MovMesh_freeBoundary_vertexNormal(X,tri,tri_bf)
%
% this function computes the unit outward normals for the boundary vertices
% of a 2D or 3D triangulation.
%
% X:    the coordinates of vertices of the mesh of size Nv-by-d.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
% tri_bf: the boundary facets for all meshes, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% V:    (output) unit outward face normals at the boundary vertices of TR,
%       of size Nv-by-d.
%       V = 0 for interior vertices.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

    [Nv,d] = size(X);
    N = size(tri, 1);
    
    if (d==1)
       error(' d = 1. MovMesh_freeBoundary_faceNormal(): d = 2 or d = 3 only.');
    end
    
    % compute outward normals for boundary segments
    V1  = MovMesh_freeBoundary_faceNormal(X,tri,tri_bf);
    
    % compute outward normals for boundary vertices
    V = zeros(Nv,d);
    for j = 1:d
    for i = 1:d
       vv = accumarray([tri_bf(:,i);Nv], [V1(:,j);0]);
       V(:,j) = V(:,j) + vv;
    end
    end
    V = bsxfun(@rdivide,V,sqrt(sum(V.^2,2)));
    V(isnan(V)) = 1/sqrt(d);

% end of MovMesh_freeBoundary_vertexNormal()
