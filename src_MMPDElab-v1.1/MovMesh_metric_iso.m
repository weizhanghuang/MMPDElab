function M = MovMesh_metric_iso(u,X,tri,tri_bf,alpha,m)
%
% usage: M = MovMesh_metric_iso(u,X,tri,tri_bf)
%        M = MovMesh_metric_iso(u,X,tri,tri_bf,alpha,m)
%
% this function computes the isotropic metric tensor M based on L2 norm or
% H1 seminorm of linear interpolation error (l = 2 and m = 0 or m = 1).
%
% u:    function values at the vertices of X, of size Nv-by-1.
% X:    coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets of the mesh, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% alpha: regularity parameter. the default value: computed through the algebraic
%       equation: int_Omega sqrt(det(M)) d x = 2 |Omega|.
%       the default value is used when alpha is not defined or empty or has
%       a nonpositive value.
% m:    m = 0 (default) for L2 norm and m = 1 for H1 seminorm.
% M:    (output) metric tensor at vertices, of size Nv-by-d*d.
%       M(i,:) is the d-by-d matrix (stored columnwise) at x_i.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   fprintf('\n--- begin MovMesh_metric_iso ---\n');

% set basic parameters

   if (nargin<6 || ~exist('m','var') || isempty(m))
      m = 0;
   end
   if ( m ~= 0 && m ~= 1)
      m = 0;
   end
   [N,d] = size(tri);
   d = d-1;
   Nv = max(max(tri));

% compute the gradient, hessian and the Frobenius norm of Hessian

   [~, Hessian] = MovMesh_GradHessianRecovery(u, X, tri, tri_bf);
   H = sqrt(sum(Hessian.^2,2))/d;
      
% compute alpha

if (nargin<5 || isempty(alpha) || alpha<=eps)

   detE = Matrix_detE(X,tri);
   detE = abs(detE)/factorial(d);
   Omega = sum(detE);
   HK = Matrix_average(H,tri).^(2*d/(d+2*(2-m)));
   alpha = (dot(detE,HK)/(eps+2*Omega)).^((d+2*(2-m))/(2*d));
   alpha = max(alpha, 0.001);
   
end
   
% compute the metric tensor 
   
   H = (1 + (1/alpha)*H).^(4/(d+2*(2-m)));
   M = repmat(reshape(eye(d),1,[]),Nv,1);
   M = bsxfun(@times,M,H);
            
   fprintf('--- end of MovMesh_metric_iso ---\n\n');

% end of MovMesh_metric_iso()
