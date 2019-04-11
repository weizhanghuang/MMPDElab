function M = MovMesh_metric(u,X,tri,tri_bf,alpha,m)
%
% usage: M = MovMesh_metric(u,X,tri,tri_bf)
%        M = MovMesh_metric(u,X,tri,tri_bf,alpha,m)
%
% this function computes the metric tensor M based on L2 norm or H1 seminorm of
% linear interpolation error (l = 2 and m = 0 or m = 1).
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

   fprintf('\n--- begin MovMesh_metric ---\n');
   
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

% compute the gradient, hessian, and |H|

   [~, Hessian] = MovMesh_GradHessianRecovery(u,X,tri,tri_bf);
   Hessian = Matrix_absA(Hessian);
      
% compute alpha using bisection

if (nargin<5 || isempty(alpha) || alpha<=eps)

   detE = Matrix_detE(X,tri);
   detE = abs(detE)/factorial(d);
   Omega = sum(detE);
   
   % estimate an upper bound (taken to be alpha for isotropic case)
   H = Matrix_average(Hessian,tri);
   HK = Matrix_trace(H).^(2*d/(d+2*(2-m)));
   b = (dot(detE,HK)/(eps+2*Omega)).^((d+2*(2-m))/(2*d));
   
   % compute alpha using bisection
   
   II = repmat(reshape(eye(d),1,[]),N,1);
   a = 0.001; 
   alpha = 0.5*(a+b); 
   while (b-a>0.001*(a+b)*0.5)
      alpha = 0.5*(a+b);
      detH = Matrix_det((1/alpha)*H + II).^((2-m)/(d+2*(2-m)));
      if (m>0)
         detH = detH.*Matrix_trace((1/alpha)*H + II).^(m*d/(d+2*(2-m)));
      end
      ff = dot(detH,detE)-2*Omega;
      if (ff<0)
         b = alpha;
      else
         a = alpha;
      end
   end
end
   
% compute the metric tensor
   
   Hessian = repmat(reshape(eye(d),1,[]),Nv,1) + (1/alpha)*Hessian;
   detHessian = Matrix_det(Hessian).^(-1/(d+2*(2-m)));
   if (m>0)
      detHessian = detHessian.*Matrix_trace(Hessian).^(m*d/(d+2*(2-m)));
   end
   M = bsxfun(@times,Hessian,detHessian);
      
   fprintf('--- end of MovMesh_metric ---\n\n');

% end of MovMesh_metric()
