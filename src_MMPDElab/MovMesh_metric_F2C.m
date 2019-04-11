function MC = MovMesh_metric_F2C(M,Tri,Tri_parent,TriC)
%
% usage: MC = MovMesh_metric_F2C(M,Tri,Tri_parent,TriC)
%
% this function computes the metric tensor on a coarse mesh
% from the metric tensor defined on a fine mesh.
%
% M:    the metric tensor defined on the fine mesh, of size Nv-by-d*d.
% Tri:  the connectivity of the fine mesh, of size N-by-(d+1).
% Tri_parent: shows which original element a new element is in,
%       of size NF-by-1. 
% TriC: the connectivity of the coarse mesh, of size NC-by-(d+1).
% MC:   (output) the metric tensor defined on the coarse mesh,
%       of size NvC-by-d*d.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   % basic parameters
   
   d = size(TriC,2)-1;
   Nv = max(max(TriC));
   N = size(TriC,1);

   % compute the elementwise averages of M

   MK = Matrix_average(M,Tri);
   
   % compute elementwise MC on the coarse mesh
   
   MCK = zeros(N,d*d);
   for j=1:d*d
      MCK(:,j) = MCK(:,j) + accumarray(Tri_parent,MK(:,j));
   end
   MCK = MCK*N/size(Tri,1);
   
   % distribute MCK to the nodes of the coarse mesh
   
   MC = zeros(Nv,d*d);
   v = zeros(Nv,1);
   for j=1:(d+1)
      for k=1:d*d
         MC(:,k) = MC(:,k) + accumarray([TriC(:,j);Nv],[MCK(:,k);0]);
      end
      v = v + accumarray([TriC(:,j);Nv], [ones(N,1);0]);
   end
   MC = bsxfun(@rdivide,MC,v);
   
% end of MovMesh_metric_F2C()
