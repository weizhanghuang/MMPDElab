function MM = MovMesh_metric_smoothing(M,ncycles,X,tri)
%
% usage: MM = MovMesh_metric_smoothing(M,ncycles,X,tri)
%
% this function smooths the metric tensor (ncycles) times by local averaging.  
%
% M:    metric tensor at vertices, of size Nv-by-d*d.
%       M(i,:) is the d-by-d matrix (stored columnwise) at x_i.
% ncycles: number of smoothing cycles.
% X:    the coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% MM:   (output) metric tensor at vertices, of size Nv-by-d*d.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   fprintf('\n--- begin MovMesh_metric_smoothing ---\n');

   [N,d] = size(tri);
   d = d-1;
   Nv = max(max(tri));
   Nm = size(M,2);
   
   MM = M;
   M1 = zeros(N,Nm,d+1);
   
for k1=1:ncycles

   wk = zeros(N,Nm);
   for j=1:(d+1)
      M1(:,:,j) = MM(tri(:,j),:);
      wk = wk + M1(:,:,j);
   end
   wk = wk/(2*d);
   for j=1:(d+1)
      M1(:,:,j) = wk + M1(:,:,j)*(d-1)/(2*d);
   end
   
   MM = zeros(Nv,Nm);
   for k=1:Nm
      for j=1:(d+1)
         v = accumarray([tri(:,j);Nv], [M1(:,k,j);0]);
         MM(:,k) = MM(:,k) + v;
      end
   end
   
   v = zeros(Nv,1);
   for j=1:(d+1)
      v = v + histc(tri(:,j),(1:Nv)');
   end
   v = 1./v;
   
   MM = bsxfun(@times,MM,v);

end
   
   fprintf('--- end of MovMesh_metric_smoothing ---\n\n');

% end of MovMesh_metric_smoothing()
