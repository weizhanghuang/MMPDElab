function [Xnew,Ih,Kmin] = MovMesh_XM(tspan,X,M,tau,tri,tri_bf,nodes_fixed, ...
                                     mmpde_int_method,dt0,abstol,Xi_ref)
%
% usage: [Xnew,Ih,Kmin] = MovMesh_XM(tspan,X,M,tau,tri,tri_bf,nodes_fixed)
%        [Xnew,Ih,Kmin] = MovMesh_XM(tspan,X,M,tau,tri,tri_bf,nodes_fixed, ...
%                                    mmpde_int_method,dt0,abstol,Xi_ref)
%
% this function computes the new physical mesh by integrating
% MMPDE (x-formulation) over tspan.
%
% tspan: a vector specifying the time interval for integration.
% X:    coordinates of vertices of the current mesh, of size Nv-by-d.
% M:    metric tensor defined at the vertices of X, of size Nv-by-d*d.
%       M(i,:) is the d-by-d matrix (stored columnwise) at x_i.
% tau:  parameter for adjusting the time scale of mesh movement.
% tri:  the connectivity for all meshes, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets for all meshes, with each row representing
%       a facet on the boundary and containing d vertex IDs,
%       of size Nbf-by-d.
% nodes_fixed: of size NvF-by-1, contains the vertex IDs which are not allowed
%       to move.
% mmpde_int_method: (optional input) method used to integrate MMPDE, chosen from
%       'ode15s' and 'ode45'. default is 'ode15s'.
% dt0:  (optional input) initial time step.
%       default: dt0 = (tspan(end)-tspan(1))/10.
% abstol: (optional input) absolute tolerance used to control time integration.
%       default: asbstol = 1e-6 for ode15s and 1e-8 for ode45.
% Xi_ref: (optional input) reference computational mesh. if not defined,
%       a mesh composed of N copies of an equilateral simplex is used.
% Xnew: (output) coordinates of vertices of the new mesh, of size Nv-by-d.
% Ih:   (optional output) the value of functional at the new mesh.
% Kmin: (optional output) the minimal element volume.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   fprintf('\n--- begin MovMesh_XM ---\n');
   
% set basic parameters
   
   if (~exist('mmpde_int_method','var') || isempty(mmpde_int_method))
      mmpde_int_method = 'ode15s';
   end
   mmpde_int_method = lower(mmpde_int_method);
   if (~strcmp(mmpde_int_method,'ode15s') && ~strcmp(mmpde_int_method,'ode45'))
      mmpde_int_method = 'ode15s';
   end  
   if (isempty(tau) || tau < eps)   
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
   if (~exist('Xi_ref','var') || isempty(Xi_ref))
      Xi_ref = [];
   end

   d = size(tri, 2) - 1;
   Nv = max(max(tri));
   
% define mesh adaptation functional

   G_fun = @MovMesh_Functional_Huang;
   %G_fun = @MovMesh_Functional_Winslow;

% for time integration

   x = reshape(X',d*Nv,1);
   xnew = x;
   reltol = 100*abstol;
   
   switch (mmpde_int_method)

   case 'ode15s'
   
      disp('    -- using ode15s');
      opts=odeset('InitialStep',dt0,'OutputFcn',@MovMesh_odeprint2, ...
                  'Jacobian',@MovMesh_jac,'AbsTol',abstol,'RelTol',reltol); 
      [~,xnew]=ode15s(@MovMesh_rhs,[tspan(1),tspan(end)],x,opts,X,M,tau, ...
                      tri,tri_bf,nodes_fixed,G_fun,Xi_ref);
      Xnew = reshape(xnew(end,:)', d, Nv)';

   case 'ode45'

      disp('    -- using ode45');
      opts=odeset('InitialStep',dt0,'OutputFcn',@MovMesh_odeprint2, ...
                  'AbsTol',abstol,'RelTol',reltol); 
      [~,xnew]=ode45(@MovMesh_rhs,[tspan(1),tspan(end)],x,opts,X,M,tau, ...
                       tri,tri_bf,nodes_fixed,G_fun,Xi_ref);
      Xnew = reshape(xnew(end,:)', d, Nv)';
   end
   
   if (nargout > 1)
      [Ih,Kmin] = MovMesh_Functional_value(Xnew,tri,tri_bf,X,M,G_fun,Xi_ref);
   end
   
   fprintf('--- end of MovMesh_XM ---\n\n');
                   
% end of MovMesh_XM()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = MovMesh_rhs(t,x,Xcurrent,M0,tau,tri,tri_bf,nodes_fixed,G_fun,Xi_ref)
%
% this function computes the rhs of the mesh equation
%
% dot(x_1, ... x_Nv) = - (1/tau) d I/ d [x_1, ... x_Nv]
%
% defined on the current mesh x.
%
% G_fun: user-supplied function defining the function G and its derivatives,
%
%       [G, GJ, GdetJ, GM] = G_fun(J, detJ, M, x, sigma).
%
%

% set basic parameters and initialization

   N = size(tri, 1);
   d = size(tri, 2) - 1;
   Nv = max(max(tri));
   Nbf = size(tri_bf,1);

   X = reshape(x,d,Nv)';
   F = zeros(Nv, d);
      
% interpolate the monitor function

   M = MovMesh_LinInterp(M0,Xcurrent,X,tri,tri_bf);
      
% computing the balancing factor b_factor

   if (strcmp(func2str(G_fun),'MovMesh_Functional_Huang')) % for Huang's functional
      p = 1.5;
      b_factor = Matrix_det(M).^((p-1)*0.5);
   else % for Winslow's functional
      b_factor = Matrix_det(M).^(1.0/d);
   end

% compute compute J, detJ, mK, xK, volK, Einv

   if (~exist('Xi_ref','var') || isempty(Xi_ref))
      switch (d)
      case 1
         Ec = 1;
      case 2
         Ec = [1, 0.5; 0, sqrt(3)*0.5];
      case 3
         Ec = [0, -2, -2; -2, 0, -2; -2, -2, 0];
      end
      Ec = (factorial(d)/abs(det(Ec)))^(1/d)*Ec; % normalize Ec
      Ec = Ec/N^(1/d);
      Ec = reshape(Ec,1,[]);
      Ec = repmat(Ec,N,1);
   else
      [Ec,~,~,~] = Matrix_edge(Xi_ref, tri);
   end
   detEc = abs(Matrix_det(Ec));
   
   [~, Einv, detE, xK] = Matrix_edge(X, tri);
   detE = abs(detE);
   volK = detE/factorial(d);
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

   [G, GJ, GdetJ, GM] = G_fun(J, detJ, mK, xK, sigma);
   
   GdetJ = GdetJ.*detJ;
   F0 = Matrix_mult(Matrix_mult(Einv,GJ),Matrix_mult(Ec,Einv));
   F0 = bsxfun(@times,Einv,G-GdetJ) - F0;
   F0 = Matrix_AT(F0);
   FF = zeros(N,d);
   for k=1:d
   for j=1:d
      FF(:,j) = FF(:,j) - F0(:,(k-1)*d+j);
   end
   end
   FF0 = [FF, F0];
   % for GM term
   for j=2:(d+1)
      tr_Mj = Matrix_trace(Matrix_mult(GM, M(tri(:,j),:)-M(tri(:,1),:)));
      EinvRowj = zeros(N,d);
      for k=1:d
         EinvRowj(:,k) = Einv(:,d*(k-1)+j-1);
      end
      EinvRowj = bsxfun(@times,EinvRowj,tr_Mj)/(d+1);
      FF0 = FF0 + repmat(EinvRowj,1,d+1);
   end
   % multiplied by |K|   
   FF0 = bsxfun(@times,FF0,volK); 
      
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
         vn = MovMesh_freeBoundary_vertexNormal(X,tri,tri_bf);
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
         nodes_b = unique(tri_bf);
         vn = MovMesh_freeBoundary_vertexNormal(X,tri,tri_bf);
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

function JAC = MovMesh_jac(t,x,Xcurrent,M0,tau,tri,tri_bf,nodes_fixed,G_fun,Xi_ref)
%
% this function computes the Jacobian matrix for the rhs of the mesh equation
%
% dot(x_1, ... x_Nv) = - (1/tau) d I/ d [x_1, ... x_Nv]
%
% defined on the current mesh x.
%
% G_fun: user-supplied function defining the function G and its derivatives,
%
%       [G, GJ, GdetJ, GM] = G_fun(J, detJ, M, x, sigma).
%

% set basic parameters and initialization

   N = size(tri, 1);
   d = size(tri, 2) - 1;
   Nv = max(max(tri));
   Nbf = size(tri_bf,1);

   X = reshape(x,d,Nv)';
   
% interpolate the monitor function

   M = MovMesh_LinInterp(M0,Xcurrent,X,tri,tri_bf);
   
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
   
   if (~exist('Xi_ref','var') || isempty(Xi_ref))
      switch (d)
      case 1
         Ec = 1;
      case 2
         Ec = [1, 0.5; 0, sqrt(3)*0.5];
      case 3
         Ec = [0, -2, -2; -2, 0, -2; -2, -2, 0];
      end
      Ec = (factorial(d)/abs(det(Ec)))^(1/d)*Ec; % normalize Ec
      Ec = Ec/N^(1/d);
      Ec = reshape(Ec,1,[]);
      Ec = repmat(Ec,N,1);
   else
      [Ec,~,~,~] = Matrix_edge(Xi_ref, tri);
   end
   detEc = abs(Matrix_det(Ec));
   
   [~, Einv, detE, xK] = Matrix_edge(X, tri);
   detE = abs(detE);
   volK = detE/factorial(d);
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
   
   b_factor = reshape(repmat(b_factor',d,1),Nv*d,1);
   
   % compute sigma

   sigma = dot(volK,Matrix_det(mK))+eps;

   % both F0 and FF are defined on elements

   [G, GJ, GdetJ, GM] = G_fun(J, detJ, mK, xK, sigma);
   
   GdetJ = GdetJ.*detJ;
   F0 = Matrix_mult(Matrix_mult(Einv,GJ),Matrix_mult(Ec,Einv));
   F0 = bsxfun(@times,Einv,G-GdetJ) - F0;
   F0 = Matrix_AT(F0);
   FF = zeros(N,d);
   for k=1:d
   for j=1:d
      FF(:,j) = FF(:,j) - F0(:,(k-1)*d+j);
   end
   end
   FF0 = [FF, F0];
   % for GM term
   for j=2:(d+1)
      tr_Mj = Matrix_trace(Matrix_mult(GM, M(tri(:,j),:)-M(tri(:,1),:)));
      EinvRowj = zeros(N,d);
      for k=1:d
         EinvRowj(:,k) = Einv(:,d*(k-1)+j-1);
      end
      EinvRowj = bsxfun(@times,EinvRowj,tr_Mj)/(d+1);
      FF0 = FF0 + repmat(EinvRowj,1,d+1);
   end
   % multiplied by |K| 
   FF0 = bsxfun(@times,FF0,volK); 
   
% compute and assemble Jacobian

ic = 1;
Indx = zeros(N,d*(d+1));
F1 = zeros(N,d*(d+1));

for j1=1:(d+1)
for k1=1:d

   [~, Einv, detE, xK, h] = Matrix_edge2(X,tri,j1,k1);
   detE = abs(detE);
   volK = detE/factorial(d);
   J = Matrix_mult(Ec, Einv);
   detJ = detEc./detE;
   
   [G, GJ, GdetJ, GM] = G_fun(J, detJ, mK, xK, sigma);
   
   GdetJ = GdetJ.*detJ;
   F1 = Matrix_mult(Matrix_mult(Einv,GJ),Matrix_mult(Ec,Einv));
   F1 = bsxfun(@times,Einv,G-GdetJ) - F1;
   F1 = Matrix_AT(F1);
   FF = zeros(N,d);
   for k=1:d
   for j=1:d
      FF(:,j) = FF(:,j) - F1(:,(k-1)*d+j);
   end
   end
   FF1 = [FF, F1];
   % for GM term
   for j=2:(d+1)
      tr_Mj = Matrix_trace(Matrix_mult(GM, M(tri(:,j),:)-M(tri(:,1),:)));
      EinvRowj = zeros(N,d);
      for k=1:d
         EinvRowj(:,k) = Einv(:,d*(k-1)+j-1);
      end
      EinvRowj = bsxfun(@times,EinvRowj,tr_Mj)/(d+1);
      FF1 = FF1 + repmat(EinvRowj,1,d+1);
   end
   % multiplied by |K|   
   FF1 = bsxfun(@times,FF1,volK);    
      
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
         vn = MovMesh_freeBoundary_vertexNormal(X,tri,tri_bf);
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
         vn = MovMesh_freeBoundary_vertexNormal(X,tri,tri_bf);
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

function [G, GJ, GdetJ, GM] = MovMesh_Functional_Huang(J, detJ, M, X, sigma)
%
% this function defines the derivatives of G wrt J and detJ for Huang's'
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

   G = theta*detM.*tr.^(d*p/2) + (1-2*theta)*d^(d*p/2)*detM.*(detJ./detM).^p;

   GdetJ = p*(1-2*theta)*d^(d*p/2)*detM.^(1-p).*detJ.^(p-1);

   GJ = Matrix_mult(Minv, Matrix_AT(J));
   GJ = d*p*theta*bsxfun(@times,GJ,detM.*tr.^(d*p/2-1));

   GM = Matrix_mult(J,Minv);
   GM = Matrix_mult(Matrix_AT(GM),GM);
   GM = - 0.5*theta*d*p*bsxfun(@times,GM,detM.*tr.^(d*p/2-1)) ...
      + 0.5*theta*bsxfun(@times,Minv,detM.*tr.^(d*p/2)) ...
      + (0.5-theta)*(1-p)*d^(d*p/2)*bsxfun(@times,Minv,detM.^(1-p).*detJ.^p);

% end of MovMesh_Functional_Huang()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [G, GJ, GdetJ, GM] = MovMesh_Functional_Winslow(J, detJ, M, X, sigma)
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

   G = Matrix_traceAMAT(J,Minv);

   GdetJ = zeros(N,1);

   GJ = 2*Matrix_mult(Minv, Matrix_AT(J));

   GM = Matrix_mult(J,Minv);
   GM = -Matrix_mult(Matrix_AT(GM),GM);

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

function [Ih,Kmin] = MovMesh_Functional_value(X, tri, tri_bf, Xcurrent, ...
                                              M0, G_fun, Xi_ref)
%
% this function computes the functional value and minimal element volume
%

% set basic parameters and initialization

   N = size(tri, 1);
   d = size(tri, 2) - 1;
   
% interpolate the monitor function

   M = MovMesh_LinInterp(M0,Xcurrent,X,tri,tri_bf);
   
% compute compute J, detJ, mK, xK, volK, Einv

   if (~exist('Xi_ref','var') || isempty(Xi_ref))
      switch (d)
      case 1
         Ec = 1;
      case 2
         Ec = [1, 0.5; 0, sqrt(3)*0.5];
      case 3
         Ec = [0, -2, -2; -2, 0, -2; -2, -2, 0];
      end
      Ec = (factorial(d)/abs(det(Ec)))^(1/d)*Ec; % normalize Ec
      Ec = Ec/N^(1/d);
      Ec = reshape(Ec,1,[]);
      Ec = repmat(Ec,N,1);
   else
      [Ec,~,~,~] = Matrix_edge(Xi_ref, tri);
   end
   detEc = abs(Matrix_det(Ec));
   
   [~, Einv, detE, xK] = Matrix_edge(X, tri);
   detE = abs(detE);
   volK = detE/factorial(d);
   J = Matrix_mult(Ec, Einv);
   detJ = detEc./detE;
   mK = Matrix_average(M,tri);
   
% compute sigma

   sigma = dot(volK,Matrix_det(mK))+eps;
   
% compute G

    [G, ~, ~, ~] = G_fun(J, detJ, mK, xK, sigma);
   
% compute I_h and Kmin

   Ih = dot(G,volK);
   Kmin = min(volK);

% end of MovMesh_Functional_value()

