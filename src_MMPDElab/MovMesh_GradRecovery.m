function Grad = MovMesh_GradRecovery(u,X,tri,tri_bf)
%
% usage: Grad = MovMesh_GradRecovery(u,X,tri,tri_bf)
%
% this function computes the gradient of function u at the vertices using volume
% averaging.
%
% u:    function values at the vertices of X, of size Nv-by-npde.
% X:    coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets of the mesh, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% Grad: (output) gradient at vertices, of size Nv-by-d*npde,
%       Grad(:,(i-1)*d+(1:d)) = \nabla u_i, i = 1, ..., npde
%       for example, d = 3: Grad = [u1_x,u1_y,u1_z,u2_x,u2_y,u2_z,...]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

% set basic parameters

   N = size(tri, 1);
   d = size(tri, 2) - 1;
   Nv = max(max(tri));
   npde = size(u,2);
   
   Grad = zeros(Nv, d*npde);
      
% compute compute volume of elements
   
   detE = Matrix_detE(X,tri);
   detE = abs(detE)/factorial(d);
   
% compute Grad at elements:

  dU = MovMesh_GradKRecovery(u,X,tri,tri_bf);
      
% compute gadient of u at nodes using volume-weighted average

   dU = bsxfun(@times,dU,detE);
   for i=1:npde
   for j=1:(d+1)
   for k=1:d
       v = accumarray([tri(:,j);Nv], [dU(:,d*(i-1)+k);0]);
       Grad(:,d*(i-1)+k) = Grad(:,d*(i-1)+k) + v;
    end
    end
    end
    v = zeros(Nv,1);
    for i=1:(d+1)
        v = v + accumarray([tri(:,i);Nv], [detE;0]);
    end
    v = 1./v;    
    Grad = bsxfun(@times,Grad,v);
         
% end of MovMesh_GradRecovery()

