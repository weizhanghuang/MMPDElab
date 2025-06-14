function M = MovMesh_metric_arclength(u,X,tri,tri_bf)
%
% usage: M = MovMesh_metric_arclength(u,X,tri,tri_bf)
%
% this function computes the arclength metric tensor.
%
% u:    function values at the vertices of X, of size Nv-by-1.
% X:    coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets of the mesh, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% M:    (output) metric tensor at vertices, of size Nv-by-d*d.
%       M(i,:) is the d-by-d matrix (stored columnwise) at x_i.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   fprintf('\n--- begin MovMesh_metric_arclength ---\n');

% set basic parameters

   [N,d] = size(tri);
   d = d-1;
   Nv = max(max(tri));

% compute the gradient

   Grad = MovMesh_GradRecovery(u,X,tri,tri_bf);
   
% compute the arclength metric tensor

   M = repmat(reshape(eye(d,d),1,[]),Nv,1);
   M = bsxfun(@times,M,sqrt(dot(Grad,Grad,2)+1));
         
   fprintf('--- end of MovMesh_metric_arclength ---\n\n');

% end of MovMesh_metric_arclength()
