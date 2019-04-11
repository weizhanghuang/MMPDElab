function [Xnew,Ih,Kmin] = MovMesh(tspan,Xi_ref,X,M,tau,tri,tri_bf,nodes_fixed, ...
                                  mmpde_int_method,dt0,abstol)
%
% usage: [Xnew,Ih,Kmin] = MovMesh(tspan,Xi_ref,X,M,tau,tri,tri_bf,nodes_fixed)
%        [Xnew,Ih,Kmin] = MovMesh(tspan,Xi_ref,X,M,tau,tri,tri_bf,nodes_fixed, ...
%                                 mmpde_int_method,dt0,abstol)
%
% this function computes the new physical mesh by integrating
% MMPDE (xi-formulation) over tspan.
%
% all meshes, Xi_ref, X should have the same N, Nv, connectivity (tri),
% boundary facets (tri_bf), and fixed nodes (nodes_fixed).
%
% tspan: a vector specifying the time interval for integration.
% Xi_ref: the coordinates of vertices of the reference computational mesh,
%       of size Nv-by-d.
% X:    the coordinates of vertices of the current mesh, of size Nv-by-d.
% M:    the metric tensor defined at the vertices of X, of size Nv-by-d*d.
%       M(i,:) is the d-by-d matrix (stored columnwise) at x_i.
% tau:  parameter for adjusting the time scale of mesh movement.
% tri:  the connectivity for all meshes, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets for all meshes, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% nodes_fixed: of size NvF-by-1, contains the vertex IDs which are not allowed
%       to move.
% mmpde_int_method: (optional input) method used to integrate MMPDE, chosen from
%       'ode15s' and 'ode45'. default is 'ode15s'.
% dt0:  (optional input) initial time step.
%       default: dt0 = (tspan(end)-tspan(1))/10.
% abstol: (optional input) absolute tolerance used to control time integration.
%       default: asbstol = 1e-6 for ode15s and 1e-8 for ode45.
% Xnew: (output) the coordinates of vertices of the new mesh, of size Nv-by-d.
% Ih:   (optional output) the value of functional at the new mesh.
% Kmin: (optional output) the minimal element volume.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   fprintf('\n--- begin MovMesh ---\n');
   
% set basic parameters
   
   if (~exist('mmpde_int_method','var') || isempty(mmpde_int_method))
      mmpde_int_method = 'ode15s';
   end
   mmpde_int_method = lower(mmpde_int_method);
   if (~strcmp(mmpde_int_method,'ode15s') && ~strcmp(mmpde_int_method,'ode45'))
      mmpde_int_method = 'ode15s';
   end
   if (~exist('tau','var') || isempty(tau) || tau < eps)   
      tau = 1e-2;   
   end
   if (tau < 0)
      tau = 1e-2;
   end
   if (~exist('dt0','var') || isempty(dt0) || dt0 < eps)
      dt0 = (tspan(end)-tspan(1))/10; 
   end
   if (dt0 < 0)
      dt0 = (tspan(end)-tspan(1))/10; 
   end
   if (~exist('abstol','var') || isempty(abstol) || abstol < eps)
      abstol = 1e-6;
      if (strcmp(mmpde_int_method,'ode45'))
         abstol = 1e-8;
      end
   end

   d = size(tri, 2) - 1;
   Nv = max(max(tri));
   
% define mesh adaptation functional

   G_fun = @MovMesh_Functional_Huang;
   %G_fun = @MovMesh_Functional_Winslow;

% for time integration

   xi = reshape(Xi_ref',d*Nv,1);
   reltol = 100*abstol;

   switch (mmpde_int_method)

   case 'ode15s'
   
      disp('    -- using ode15s');
      opts=odeset('InitialStep',dt0,'OutputFcn',@MovMesh_odeprint2, ...
                  'Jacobian',@MovMesh_jac,'AbsTol',abstol,'RelTol',reltol); 
      [~,xinew]=ode15s(@MovMesh_rhs,[tspan(1),tspan(end)],xi,opts,X, ...
                       M,tau,tri,tri_bf,nodes_fixed,G_fun);       
      Xinew = reshape(xinew(end,:)', d, Nv)';

   case 'ode45'

      disp('    -- using ode45');
      opts=odeset('InitialStep',dt0,'OutputFcn',@MovMesh_odeprint2, ...
                  'AbsTol',abstol,'RelTol',reltol);
      [~,xinew]=ode45(@MovMesh_rhs,[tspan(1),tspan(end)],xi,opts,X,M,tau, ...
                       tri,tri_bf,nodes_fixed,G_fun);
      Xinew = reshape(xinew(end,:)', d, Nv)';
   end
             
% compute the physical mesh using linear interpolation

   if (d == 1)
      Xnew = MovMesh_LinInterp(X,Xinew,Xi_ref,tri,tri_bf,true);
   else
      % compute the volume of the convex hull
      DT = delaunayTriangulation(Xinew);
      [~,vol] = convexHull(DT);
      [~,~,detE,~] = Matrix_edge(Xinew,tri);
      vol2 = sum(abs(detE))/factorial(d);
      if (abs(vol)<=vol2+100*eps) % convex, using delauney class
         Xnew = MovMesh_LinInterp(X,Xinew,Xi_ref,tri,tri_bf,true);
      else % not convex, use triangulation class
         Xnew = MovMesh_LinInterp(X,Xinew,Xi_ref,tri,tri_bf,false);
         % re-compute the boundary nodes
         TR = triangulation(tri,Xinew);
         % find the elements containing the boundary facets
         if (d == 2)
            bfElements = edgeAttachments(TR,tri_bf);
         else % d = 3
            Elem1 = edgeAttachments(TR,tri_bf(:,1),tri_bf(:,2));
            bfElements = edgeAttachments(TR,tri_bf(:,2),tri_bf(:,3));
            bfElements = cellfun(@intersect,Elem1,bfElements,'UniformOutput',false);
            Elem1 = edgeAttachments(TR,tri_bf(:,1),tri_bf(:,3));
            bfElements = cellfun(@intersect,Elem1,bfElements,'UniformOutput',false);
         end
         bfE = cell2mat(bfElements);
         n_bfE = size(bfE,1);
         B = zeros(n_bfE,d+1);
         nodes = setdiff(unique(tri_bf),nodes_fixed);
         n_nodes = size(nodes,1);
         for i=1:n_nodes
            node = nodes(i);
            B = cartesianToBarycentric(TR,bfE,repmat(Xi_ref(node,:),n_bfE,1));
            [~,ii] = min(abs(sum(B,2)-sum(abs(B),2)));
            Xnew(node,:) = zeros(1,d);
            for j=1:d+1
              Xnew(node,:) = Xnew(node,:) + B(ii,j)*X(tri(bfE(ii),j),:);
            end
         end
      end
   end
      
   % for optinal output
   
   if (nargout > 1)
      [Ih,Kmin] = MovMesh_Functional_value(Xinew,tri,X,M,G_fun);
   end
   
   fprintf('--- end of MovMesh ---\n\n');
                   
% end of MovMesh()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = MovMesh_rhs(t,xi,X,M,tau,tri,tri_bf,nodes_fixed,G_fun)
%
% this function computes the rhs of the mesh equations
%
% dot(xi_1, ... xi_Nv) = - (1/tau) d I/ d [xi_1, ... xi_Nv]
%
% defined on the current meshes (xi, X).
%
% G_fun: user-supplied function defining the function G and its derivatives,
%
%       [G, GJ, GdetJ] = G_fun(J, detJ, M, x, sigma)
%
%

% set basic parameters and initialization

   N = size(tri, 1);
   d = size(tri, 2) - 1;
   Nv = max(max(tri));
   Nbf = size(tri_bf,1);

   XI = reshape(xi,d,Nv)';
   F = zeros(Nv, d);
   
% computing the balancing factor b_factor

   if (strcmp(func2str(G_fun),'MovMesh_Functional_Huang')) % for Huang's functional
      p = 1.5;
      b_factor = Matrix_det(M).^((p-1)*0.5);
   else % for Winslow's functional
      b_factor = Matrix_det(M).^(1.0/d);
   end

% compute compute J, detJ, mK, xK, volK, Einv, Ecinv
   
   [~, Einv, detE, xK] = Matrix_edge(X, tri);
   detE = abs(detE);
   volK = detE/factorial(d);
   [Ec, Ecinv, detEc, ~] = Matrix_edge(XI, tri);
   detEc = abs(detEc);
   J = Matrix_mult(Ec, Einv);
   detJ = detEc./detE;
   mK = Matrix_average(M,tri);
   
% multiply b_factor with volume-based factor, i.e., tau ==> tau \omega_i

%{
b_factor_volume = false;
if (b_factor_volume) 
   omegaj = zeros(Nv,1);
   for j=1:(d+1)
      v = accumarray([tri(:,j);Nv], [volK;0]);
      omegaj = omegaj + v; 
   end
   b_factor = b_factor./omegaj;
end
%}
   
% compute sigma

   sigma = dot(volK,Matrix_det(mK))+eps;

% compute rhs and assemble them

   % both F0 and FF are defined on elements

   [~, GJ, GdetJ] = G_fun(J, detJ, mK, xK, sigma);
   GdetJ = GdetJ.*detJ;
   F0 = Matrix_mult(Einv,GJ) + bsxfun(@times,Ecinv,GdetJ);
   F0 = bsxfun(@times,F0,volK);   
   F0 = Matrix_AT(F0);
   FF = zeros(N,d);
   for k=1:d
   for j=1:d
      FF(:,j) = FF(:,j) - F0(:,(k-1)*d+j);
   end
   end
   FF0 = [FF, F0];
      
   % F is defined on nodes
  
   for j=1:(d+1)
   for k=1:d
      v = accumarray([tri(:,j);Nv], [FF0(:,d*(j-1)+k);0]);
      F(:,k) = F(:,k) + v; 
   end
   end

   F = -(1/tau)*bsxfun(@times,F,b_factor);
      
% modify forces at boundary points
   
   switch (d)
   
   case 1
   
      F(tri_bf,:) = 0;
   
   case 2
   
      fixed = false;
      if (fixed)
         F(tri_bf,:) = 0;
      else
         % project the velocities onto the boundary
         vn = MovMesh_freeBoundary_vertexNormal(XI,tri,tri_bf);
         nodes_b = unique(tri_bf);
         w = dot(vn(nodes_b,:), F(nodes_b,:),2);
         F(nodes_b,:) = F(nodes_b,:)-bsxfun(@times,vn(nodes_b,:),w);
      end
      
   case 3
   
      fixed = false;
      if (fixed) 
         F(tri_bf,:) = 0;
      else
         % project the velocities onto the boundary
         vn = MovMesh_freeBoundary_vertexNormal(XI,tri,tri_bf);
         nodes_b = unique(tri_bf);
         w = dot(vn(nodes_b,:), F(nodes_b,:),2);
         F(nodes_b,:) = F(nodes_b,:)-bsxfun(@times,vn(nodes_b,:),w);
      end      
   end

% modify rhs for fixed nodes

   nvf = length(nodes_fixed);
   F(nodes_fixed,:) = zeros(nvf,d);
   
% convert F to vector

   F = reshape(F',d*Nv,1);

% end of MovMesh_rhs()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function JAC = MovMesh_jac(t,xi,X,M,tau,tri,tri_bf,nodes_fixed,G_fun)
%
% this function computes the Jacobian matrix for the rhs of the mesh equations
%
% dot(xi_1, ... xi_Nv) = - (1/tau) d I/ d [xi_1, ... xi_Nv]
%
% defined on the current meshes (xi, X).
%
% G_fun: user-supplied function defining the function G and its derivatives,
%
%       [G, GJ, GdetJ] = G_fun(J, detJ, M, x, sigma)
%

% set basic parameter and initialize

   N = size(tri, 1);
   d = size(tri, 2) - 1;
   Nv = max(max(tri));
   Nbf = size(tri_bf,1);

   XI = reshape(xi,d,Nv)';

   nA = (d*(d+1))^2;
   II = zeros(N*nA,1);
   JJ = zeros(N*nA,1);
   JAC_nz = zeros(N*nA,1);

% computing the balancing factor b_factor
     
   if (strcmp(func2str(G_fun),'MovMesh_Functional_Huang')) % for Huang's functional
      p = 1.5;
      b_factor = Matrix_det(M).^((p-1)*0.5);
   else % for Winslow's functional
      b_factor = Matrix_det(M).^(1.0/d);
   end

% compute rhs

   % compute compute J, detJ, mK, xK, volK, Einv, Ecinv   
   [~, Einv, detE, xK] = Matrix_edge(X, tri);
   detE = abs(detE);
   volK = detE/factorial(d);
   [Ec, Ecinv, detEc, ~] = Matrix_edge(XI, tri);
   detEc = abs(detEc);
   J = Matrix_mult(Ec, Einv);
   detJ = detEc./detE;
   mK = Matrix_average(M, tri);
   
% multiply b_factor with volume-based factor, i.e., tau ==> tau \omega_i

%{
b_factor_volume = false;
if (b_factor_volume) 
   omegaj = zeros(Nv,1);
   for j=1:(d+1)
      v = accumarray([tri(:,j);Nv], [volK;0]);
      omegaj = omegaj + v; 
   end
   b_factor = b_factor./omegaj;
end
%}
   
   b_factor = reshape(repmat(b_factor',d,1),Nv*d,1);
   
   % compute sigma

   sigma = dot(volK,Matrix_det(mK))+eps;
   
   % both F0 and FF are defined on elements
   
   [~, GJ, GdetJ] = G_fun(J, detJ, mK, xK, sigma);
   GdetJ = GdetJ.*detJ;
   F0 = Matrix_mult(Einv,GJ) + bsxfun(@times,Ecinv,GdetJ);   
   F0 = bsxfun(@times,F0,volK);
   F0 = Matrix_AT(F0);
   FF = zeros(N,d);
   for k=1:d
   for j=1:d
      FF(:,j) = FF(:,j) - F0(:,(k-1)*d+j);
   end
   end
   FF0 = [FF, F0];
   
% compute and assemble Jacobian

ic = 1;
Indx = zeros(N,d*(d+1));
F1 = zeros(N,d*(d+1));

for j1=1:(d+1)
for k1=1:d

  [Ec,Ecinv,detEc,~,h] = Matrix_edge2(XI,tri,j1,k1);
  detEc = abs(detEc);
  J = Matrix_mult(Ec, Einv);
  detJ = detEc./detE;
  [~, GJ, GdetJ] = G_fun(J, detJ, mK, xK, sigma);
   GdetJ = GdetJ.*detJ;
   F1 = Matrix_mult(Einv,GJ) + bsxfun(@times,Ecinv,GdetJ);
   F1 = bsxfun(@times,F1,volK);
   F1 = Matrix_AT(F1);
   FF = zeros(N,d);
   for k=1:d
   for j=1:d
      FF(:,j) = FF(:,j) - F1(:,(k-1)*d+j);
   end
   end
   FF1 = [FF, F1];
      
   h = 1./h;
   FF1 = FF1-FF0;
   FF1 = bsxfun(@times,FF1,h);
   JAC_nz(ic:(ic+N*d*(d+1)-1)) = reshape(FF1, N*d*(d+1), 1);
   for j=1:(d*(d+1))
      Indx(:,j) = d*(tri(:,j1)-1)+k1;
   end
   JJ(ic:(ic+N*d*(d+1)-1)) = reshape(Indx, N*d*(d+1), 1);
   for j=1:(d+1)
   for k=1:d
      Indx(:,(j-1)*d+k) = d*(tri(:,j)-1)+k;
   end
   end
   II(ic:(ic+N*d*(d+1)-1)) = reshape(Indx, N*d*(d+1), 1);
   ic = ic + N*d*(d+1);
end
end

% save to sparse matrix

   JAC_nz = JAC_nz.*b_factor(II);
   JAC = sparse(II, JJ, JAC_nz, Nv*d, Nv*d);
   JAC = - JAC/tau;

% use transpose of JAC since it is more efficient
% to deal with columns than rows

   JAC = JAC';
   
% modify jacobian matrix at boundary points
   
   % for the fixed node situation

   switch (d)
   
   case 1
   
      JAC(:,tri_bf) = 0;
   
   case 2
   
      fixed = false;
      if (fixed) 
         JAC(:,d*(tri_bf-1)+1) = 0;
         JAC(:,d*(tri_bf-1)+2) = 0;
      else 
         % project the velocities onto the boundary
         vn = MovMesh_freeBoundary_vertexNormal(XI,tri,tri_bf);
         nodes_b = unique(tri_bf);
         JAC1 = JAC(:,d*(nodes_b-1)+1);
         JAC2 = JAC(:,d*(nodes_b-1)+2);
         for i=1:d
            JAC(:,d*(nodes_b-1)+i) = JAC(:,d*(nodes_b-1)+i) ...
                    - bsxfun(@times,JAC1,(vn(nodes_b,i).*vn(nodes_b,1)).') ...
                    - bsxfun(@times,JAC2,(vn(nodes_b,i).*vn(nodes_b,2)).');
         end
      end
   
   case 3
      
      fixed = false;
      if (fixed)
         JAC(:,d*(tri_bf-1)+1) = 0;
         JAC(:,d*(tri_bf-1)+2) = 0;
         JAC(:,d*(tri_bf-1)+3) = 0;
      else
         % project the velocities onto the boundary
         vn = MovMesh_freeBoundary_vertexNormal(XI,tri,tri_bf);
         nodes_b = unique(tri_bf);
         JAC1 = JAC(:,d*(nodes_b-1)+1);
         JAC2 = JAC(:,d*(nodes_b-1)+2);
         JAC3 = JAC(:,d*(nodes_b-1)+3);
         for i=1:d
            JAC(:,d*(nodes_b-1)+i) = JAC(:,d*(nodes_b-1)+i) ...
                    - bsxfun(@times,JAC1,(vn(nodes_b,i).*vn(nodes_b,1)).') ...
                    - bsxfun(@times,JAC2,(vn(nodes_b,i).*vn(nodes_b,2)).') ...
                    - bsxfun(@times,JAC3,(vn(nodes_b,i).*vn(nodes_b,3)).');
         end
      end      
   end
   
% modify jacobian matrix for fixed nodes

   switch (d)
   case 1
      JAC(:,nodes_fixed) = 0;
   case 2
      JAC(:,d*(nodes_fixed-1)+1) = 0;
      JAC(:,d*(nodes_fixed-1)+2) = 0;
   case 3
      JAC(:,d*(nodes_fixed-1)+1) = 0;
      JAC(:,d*(nodes_fixed-1)+2) = 0;
      JAC(:,d*(nodes_fixed-1)+3) = 0;
   end
   
   JAC = JAC'; 
   
% end of MovMesh_jac()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G, GJ, GdetJ] = MovMesh_Functional_Huang(J, detJ, M, X, sigma)
%
% this function defines the derivatives of G wrt J and detJ for Huang's
% functional based on equidistribution and alignment conditions.
%
% J:    Jacobian matrix of xi = xi(x), of size N-by-d*d.
% detJ: the Jacobian determinant, of size N-by-1.
% M:    metric tensor, of size N-by-d*d.
% X:    physical coordinate, of size N-by-d.
% G:    (output) integrand, of size N-by-1.
% GJ:   (output) derivative of G wrt J, of size N-by-d*d.
% GdetJ: (output) derivative of G wrt detJ, of size N-by-1.
%

   d = size(X,2);

   theta = 1.0/3.0; % 0 < theta <= 0.5
   p = 1.5; % p > 1

   detM = Matrix_det(M);
   Minv = Matrix_inv(M,detM);
   tr = Matrix_traceAMAT(J,Minv);

   detM = sqrt(detM);

   G = theta*detM.*tr.^(d*p/2)+(1-2*theta)*d^(d*p/2)*detM.^(1-p).*detJ.^p;

   GdetJ = p*(1-2*theta)*d^(d*p/2)*detM.^(1-p).*detJ.^(p-1);

   GJ = Matrix_mult(Minv, Matrix_AT(J));
   GJ = d*p*theta*bsxfun(@times,GJ,detM.*tr.^(d*p/2-1));

% end of MovMesh_Functional_Huang()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G, GJ, GdetJ] = MovMesh_Functional_Winslow(J, detJ, M, X, sigma)
%
% this function defines the derivatives of G wrt J and detJ
% for the generalization of Winslow's functional of variable diffusion.
%
% J:    Jacobian matrix of xi = xi(x), of size N-by-d*d.
% detJ: the Jacobian determinant, of size N-by-1.
% M:    metric tensor, of size N-by-d*d.
% X:    physical coordinate, of size N-by-d.
% G:    (output) integrand, of size N-by-1.
% GJ:   (output) derivative of G wrt J, of size N-by-d*d.
% GdetJ: (output) derivative of G wrt detJ, of size N-by-1.
%

   N = size(J,1);

   detM = Matrix_det(M);
   Minv = Matrix_inv(M,detM);

   G = 0.5*Matrix_traceAMAT(J,Minv);

   GdetJ = zeros(N,1);
   GJ = Matrix_mult(Minv, Matrix_AT(J));

% end of MovMesh_Functional_Winslow()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function status = MovMesh_odeprint2(t,y,flag,varargin)

if nargin < 3 || isempty(flag) % odeprint(t,y) [v5 syntax] or odeprint(t,y,'')
  
    fprintf('      t = %e y = %e\n', t(end), y(1));
    
else
  switch(flag)
  case 'init'               % odeprint(tspan,y0,'init')

    %fprintf('      t = %e y = %e\n', t(end), y(1));
    
  case 'done'               % odeprint([],[],'done')
    
  end
end

status = 0;

% end of MovMesh_odeprint2()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ih,Kmin] = MovMesh_Functional_value(XI, tri, X, M, G_fun)
%
% this function computes the functional value and minimal element volume
%

% set basic parameters and initialization

   d = size(XI, 2);

% compute compute J, detJ, mK, xK, volK, Einv, Ecinv
   
   [~, Einv, detE, xK] = Matrix_edge(X, tri);
   detE = abs(detE);
   volK = detE/factorial(d);
   [Ec, ~, detEc, ~] = Matrix_edge(XI, tri);
   detEc = abs(detEc);
   J = Matrix_mult(Ec, Einv);
   detJ = detEc./detE;
   mK = Matrix_average(M,tri);
   
% compute sigma

   sigma = dot(volK,Matrix_det(mK))+eps;

% compute G

   [G, ~, ~] = G_fun(J, detJ, mK, xK, sigma);
   
% compute I_h and Kmin

   Ih = dot(G,volK);
   Kmin = min(volK);

% end of MovMesh_Functional_value()
