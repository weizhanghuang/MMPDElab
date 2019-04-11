function ex1d_heat()
%
% example for heat equation on Omega = (0,1), t in (0,1]
%
% Dirichlet BC on x = 0 & x = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 61;

   epsilon = 1e-3;
   npde = 1;
   moving_mesh = true;
   mmpde_tau = 1e-3;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   t = 0;
   tf = 1;
   dt0 = 1e-2;
   dtmax = 0.1;
   
% set the initial meshes, find the indices of the corner points and fix them

   X = linspace(0,1,jmax)';
   tri = [(1:jmax-1)',(2:jmax)'];
   tri_bf = [1;jmax];
   Nbf = size(tri_bf,1);

   Nv = jmax;
   d = 1;
   N = jmax-1;
   Xi_ref = X;

   % define nodes_fixed: corners are fixed
   nodes_fixed = [1;jmax];
   
% set initial conditions and compute the initial adjusted mesh

   % set the initial solution
   U = uexact(t, X);
   
   figure(1)
   clf
   Ue = uexact(t,X);
   plot(X,U(:,1),'-bo',X,Ue(:,1),'-r');
   axis([0 1 -1.1 1.1]);
   drawnow;

% define PDE system and BCs

   % all bcs are dirichlet so no need for marking boundary segments
   pdedef.bfMark = ones(Nbf,1);
   pdedef.bftype = ones(Nbf,npde);

   pdedef.volumeInt = @pdedef_volumeInt;
   pdedef.boundaryInt = @pdedef_boundaryInt;
   pdedef.dirichletRes = @pdedef_dirichletRes;

% perform integration (MP)

   DT = zeros(20000,4);
   tcpu = cputime;
   err_total = 0.0;
   n = 0;
   dt = dt0;
   
   while true
      
      % move the mesh
      
      if (moving_mesh)
         M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
         M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
         Xnew = MovMesh([t,t+dt],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
      else
         Xnew = X;
      end

      Xdot = (Xnew-X)/dt;
      
      % integrate physical PDEs
      
      [Unew,dt0,dt1] = MovFEM(t,dt,U,X,Xdot,tri,tri_bf,pdedef);
      
      % update

      X = X + dt0*Xdot;
      U = Unew;
      n = n + 1;
      
      t = t + dt0;
      dt = min(dtmax,dt1);
      if (t+dt>tf), dt=tf-t; end

      err = MovFEM_Error_P1L2(@uexact,t,X,U,tri,tri_bf);
      err_total=err_total+dt0*err;
      fprintf('--- n = %d  t = %e dt0 = %e dt1 = %e error = %e\n', ...
                     n,t,dt0,dt1,err_total);
      
      DT(n,:) = [t, dt0, err_total, err]; 
      
      figure(2)
      Ue = uexact(t,X);
      plot(X,U(:,1),'-bo',X,Ue(:,1),'-r');
      axis([0 1 -1.1 1.1]);
      drawnow;
      
      if (t>=tf-100*eps || dt < 100*eps), break; end
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
      
   figure(4)
   clf
   semilogy(DT(:,1),DT(:,2),'-o');
   title('dt')
   
   figure(5)
   clf
   semilogy(DT(:,1),DT(:,3),'-x',DT(:,1),DT(:,4),'-');
   title('err_total, err')
   
   Ue = uexact(t,X);
   fprintf('(Nv, N, max err, L2 err) = %d %d %e %e\n', ...
                              Nv, N, norm(Ue-U,Inf), err_total);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uexact(t, x)

   u(:,1) = cos(2*pi*x(:,1))*exp(-t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)

    F = - cos(2*pi*x(:,1))*exp(-t) + 4*pi*pi*cos(2*pi*x(:,1))*exp(-t);
    F = ut(:,1).*v(:) + du(:,1).*dv(:,1) - F.*v(:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   G = zeros(size(x,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u(:,1)-uexact(t,x);
