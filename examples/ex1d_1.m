function ex1d_1()
%
% example for adaptive mesh generation for given functions for intervals
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   % mmpde_int_method: 'ode15s' (default) or 'ode45'
   % mmpde_dt0: (tspan(end-1)-tspan(1))/5 or other positive numbers
   % mmpde_abstol: 1e-6 (default for ode15s), 1e-8 (default for ode15s)
   % method: 'MovMesh', 'MovMesh_XM', or 'MovMesh_X'

   mmpde_tau = 1e-2;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   mmpde_int_method = [];
   mmpde_dt0 = [];
   mmpde_tf = 0.1;
   mmpde_abstol = [];
   method = 'MovMesh';
   nn = 10;

% set the initial meshes, find the indices of the corner points and fix them 

   jmax = 61;  

   X = linspace(0,1,jmax)';
   tri = [(1:jmax-1)',(2:jmax)'];
   tri_bf = [1;jmax];
   
   Nv = jmax;
   d = 1;
   N = jmax-1;
   Xi_ref = X;

   % define nodes_fixed: corners are fixed
   nodes_fixed = [1;jmax];
   
   % perturb the mesh
%{      
   h = (1/N)^(1/d);
   dX = (2*rand(Nv,d)-1)*0.5*h*0.8;
   dX(tri_bf,:) = 0;
   X = X + dX;
%}
   
   figure(1)
   clf
   U = uexact(0, X);
   plot(X,U,'r-o');
   axis([-0.1 1.1 -1.4 1.1]);
   hold on
   plot(X,ones(Nv,1)*(-1.2),'+r');
   hold off
   drawnow;
   
% perform integration (to generate the adaptive mesh)

   TT = zeros(nn,1);
   Ih = zeros(nn,1);
   Kmin = zeros(nn,1);
   
   tcpu = cputime;

   for n=1:nn

      % compute metric tensor      
      U = uexact(0, X);
      M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
      M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);

      if (n==1)
         [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
      end
      
      switch (method)
      case 'MovMesh' % xi-formulation
         [Xnew,Ih(n),Kmin(n)] = MovMesh([0,mmpde_tf],Xi_ref,X,M,mmpde_tau,tri,tri_bf, ...
                        nodes_fixed,mmpde_int_method,mmpde_dt0,mmpde_abstol);
      case 'MovMesh_XM' % x-formulation
         [Xnew,Ih(n),Kmin(n)] = MovMesh_XM([0,mmpde_tf],X,M,mmpde_tau,tri,tri_bf, ...
                       nodes_fixed,mmpde_int_method,mmpde_dt0,mmpde_abstol);
      case 'MovMesh_X' % x-formulation with M = I (no adaptation)
         [Xnew,Ih(n),Kmin(n)] = MovMesh_X([0,mmpde_tf],X,mmpde_tau,tri,tri_bf, ...
                       nodes_fixed,mmpde_int_method,mmpde_dt0,mmpde_abstol);
      end

      fprintf('--- n = %d  Ih = %e\n', n, Ih(n));
      
      TT(n) = n;
      X = Xnew;
      
      figure(2)
      U = uexact(0, X);
      plot(X,U,'-o');
      axis([-0.1 1.1 -1.4 1.1]);
      hold on
      plot(X,ones(Nv,1)*(-1.2),'+r');
      hold off
      drawnow;
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);

% output

   U = uexact(0,X);
   err = MovFEM_Error_P1L2(@uexact,0,X,U,tri,tri_bf);
   fprintf('(Nv, N, error) = %d %d %e\n', Nv, N, err);
   fprintf('initial mesh (Qgeo, Qeq, Qali) = %e %e %e\n', Qgeo, Qeq, Qali)
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('final   mesh (Qgeo, Qeq, Qali) = %e %e %e\n', Qgeo, Qeq, Qali)

   figure(4)
   clf
   plot(TT,Ih,'r-o');
   title('Ih')

   figure(5)
   clf
   semilogy(TT(1:end-1),abs(Ih(2:end)-Ih(1:end-1)),'-o');
   title('Ih - diff')

   figure(6)
   clf
   semilogy(TT,Kmin,'r-o');
   title('Kmin')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uexact(t, x)

u = zeros(size(x,1),1);

example = 2;

switch (example)

case 1

   u = sin(2*pi.*x(:,1));

case 2

   R = 30;
   u = tanh(R*(x(:,1)-0.5));
   
case 3

   R = 30;
   u = tanh(R*x);
      
end

