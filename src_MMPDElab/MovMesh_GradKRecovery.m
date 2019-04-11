function Grad = MovMesh_GradKRecovery(u,X,tri,tri_bf)
%
% usage: Grad = MovMesh_GradKRecovery(u,X,tri,tri_bf)
%
% this function computes the gradient of function u on elements. 
%
% u:    function values at the vertices of X, of size Nv-by-npde.
% X:    coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  connectivity of the mesh, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets of the mesh, with each row representing
%       a facet on the boundary and containing d vertex IDs,
%       of size Nbf-by-d.
% Grad: (output) gradient on elements, of size N-by-d*npde,
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
      
% compute compute Einv and Phi
   
   [~, Einv, ~, ~] = Matrix_edge(X, tri);
   Phi = zeros(N,d*(d+1));
   for j=1:d
   for k=1:d
      Phi(:,j) = Phi(:,j) - Einv(:,(j-1)*d+k);
   end
   end
   for j=1:d
   for k=1:d
      Phi(:,d+(j-1)*d+k) = Einv(:,(k-1)*d+j);
   end
   end
      
% compute Grad at elements: Grad(:,(i-1)*d+(1:d)) = \nabla u_i

    Grad = zeros(N,d*npde);
    for i=1:npde
    for j=1:(d+1)
    for k=1:d
        Grad(:,d*(i-1)+k)=Grad(:,d*(i-1)+k)+Phi(:,(j-1)*d+k).*u(tri(:,j),i);
    end
    end
    end

% end of MovMesh_GradKRecovery()
