function V = MovMesh_freeBoundary_faceNormal(X,tri,tri_bf)
%
% usage: V = MovMesh_freeBoundary_faceNormal(X,tri,tri_bf)
%
% this function computes the unit outward normals for the boundary facets
% of a 2D or 3D triangulation.
%
% X:    the coordinates of vertices of the mesh of size Nv-by-d.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
% tri_bf: the boundary facets for all meshes, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% V:    (output) unit outward normals of the boundary facets, of size Nbf-by-d.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

    [Nv,d] = size(X);
    N = size(tri, 1);
    
    switch (d)
    case 1
       error(' d = 1. MovMesh_freeBoundary_faceNormal(): d = 2 or d = 3 only.');
    case 2
       % compute normals for boundary segments
       V = X(tri_bf(:,1),:)-X(tri_bf(:,2),:);
       V = bsxfun(@rdivide,V,sqrt(sum(V.^2,2)));
       V(isnan(V)) = 1/sqrt(d);
       V = [-V(:,2),V(:,1)];
       
       % determine outward normals
       TR = triangulation(tri,X);
       ID = edgeAttachments(TR,tri_bf);
       C = incenter(TR,cell2mat(ID));
       C = 0.5*(X(tri_bf(:,1),:)+X(tri_bf(:,2),:))-C;
       C(:,1) = V(:,1).*C(:,1)+V(:,2).*C(:,2);
       V(C(:,1)<0,:) = - V(C(:,1)<0,:);
    case 3
       V1 = X(tri_bf(:,1),:)-X(tri_bf(:,2),:);
       V2 = X(tri_bf(:,1),:)-X(tri_bf(:,3),:);
       V = [V1(:,2).*V2(:,3)-V1(:,3).*V2(:,2), ...
           -V1(:,1).*V2(:,3)+V1(:,3).*V2(:,1), ...
            V1(:,1).*V2(:,2)-V1(:,2).*V2(:,1)];
       V = bsxfun(@rdivide,V,sqrt(sum(V.^2,2)));
       V(isnan(V)) = 1/sqrt(d);
       
       % determine outward normals
       TR = triangulation(tri,X);
       Elem1 = edgeAttachments(TR,tri_bf(:,1),tri_bf(:,2));
       Elem2 = edgeAttachments(TR,tri_bf(:,2),tri_bf(:,3));
       Elem = cellfun(@intersect,Elem1,Elem2,'UniformOutput',false);
       Elem1 = edgeAttachments(TR,tri_bf(:,1),tri_bf(:,3));
       Elem2 = cellfun(@intersect,Elem1,Elem,'UniformOutput',false);
       C = incenter(TR,cell2mat(Elem2));
       C = (X(tri_bf(:,1),:)+X(tri_bf(:,2),:)+X(tri_bf(:,3),:))/3-C;
       C(:,1) = V(:,1).*C(:,1)+V(:,2).*C(:,2)+V(:,3).*C(:,3);
       V(C(:,1)<0,:) = - V(C(:,1)<0,:);
    end

% end of MovMesh_freeBoundary_faceNormal()
