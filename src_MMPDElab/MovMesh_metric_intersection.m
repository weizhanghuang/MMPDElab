function M = MovMesh_metric_intersection(M1,M2)
%
% usage: M = MovMesh_metric_intersection(M1,M2)
%
% this function computes the intersection of two nonsingular, symmetric
% matrices.
%
% M1 and M2: metric tensors at vertices, of size Nv-by-d*d.
% M:    (output) metric tensor at vertices, of size Nv-by-d*d.
%       M(i,:) is the d-by-d matrix.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   [Nv,d] = size(M1);
   d = sqrt(d);
   M = zeros(Nv,d*d);
   
   for i=1:Nv
      U = chol(reshape(M1(i,:),d,d));
      Uinv = inv(U);
      MK = Uinv'*reshape(M2(i,:),d,d)*Uinv;
      [V,D] = eig(MK);
      D(1:(d+1):end) = max(D(1:(d+1):end),1);
      U = U'*V;
      MK = U*D*U';
      M(i,:) = reshape(MK,1,[]);
   end

% end of MovMesh_metric_intersection()
