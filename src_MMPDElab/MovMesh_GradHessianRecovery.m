function [Grad,Hessian] = MovMesh_GradHessianRecovery(u,X,tri,tri_bf)
%
% usage: [Grad,Hessian] = MovMesh_GradHessianRecovery(u,X,tri,tri_bf)
%
% this function computes the gradient and Hessian of function u at
% the vertices using centroid-vortex-centroid-vertex volume-weighted average.
%
% u:    function values at the vertices of X, of size Nv-by-1.
% X:    coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  connectivity of the mesh, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets of the mesh, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% Grad: (output) gradient at vertices, of size Nv-by-d,
%       Grad(i,:) = [Grad_x Grad_y Grad_z](x_i).
% Hessian: (output) Hessian at vertices, of size Nv-by-d*d,
%       Hessian(i,:) is the d-by-d matrix (stored columnwise) at x_i.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)
   
% set basic parameters

   disp('    -- using CVCV-average');

   d = size(tri, 2) - 1;
   Nv = max(max(tri));
   
   Grad = MovMesh_GradRecovery(u,X,tri,tri_bf);
   Hessian = zeros(Nv, d*d);
   for i=1:d
      Hessian(:,d*(i-1)+(1:d)) = MovMesh_GradRecovery(Grad(:,i),X,tri,tri_bf);
   end
   Hessian = 0.5*(Hessian+Matrix_AT(Hessian));
   Hessian(tri_bf,:) = 0;  % throw away hessian at boundary points
         
% end of MovMesh_GradHessianRecovery()

