function ex2d_ut()
%
% example for u_t = 0,  Omega = (0,1) x (0,1), t in (0,1]
% 
% no boundary conditions ("natural boundary conditions")
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 31;

   npde = 1;
   moving_mesh = true;
   mmpde_tau = 1e-2;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   t = 0;
   tf = 1;
   dt0 = 1e-2;
   dtmax = 0.1;
   
% set the initial meshes, find the indices of the corner points and fix them

   kmax = jmax;
   [X,tri] = MovMesh_rect2tri(linspace(0,1,jmax),linspace(0,1,kmax),1);
   TR = triangulation(tri,X);
   tri_bf = freeBoundary(TR);
   Nbf = length(tri_bf);

   [Nv,d] = size(X);
   N = size(tri, 1);
   Xi_ref = X;

   % find the indices of the corner points and fix them
   corners = [0, 0; 1, 0; 1, 1; 0, 1];
   [~,nodes_fixed] = ismembertol(corners,Xi_ref,1e-10,'ByRows',true);
   % nodes_fixed = unique(tri_bf); % for fixing all boundary nodes

% set initial conditions

   % set the initial solution
   U = uexact(t,X);
      
   figure(1)
   clf
   triplot(tri,X(:,1),X(:,2),'Color','r')
   axis([0 1 0 1]);
   axis square;

% define PDE system and BCs
   
   % no need to mark the boundary segments since neumann bcs for all
   % boundary segments.
   
   pdedef.bfMark = ones(Nbf,1); 
   pdedef.bftype = zeros(Nbf,npde);

   pdedef.volumeInt = @pdedef_volumeInt;
   pdedef.boundaryInt = @pdedef_boundaryInt;
   pdedef.dirichletRes = @pdedef_dirichletRes;
   
% perform integration (MP)

   dt = dt0;
   DT = zeros(20000,2);
   tcpu = cputime;
   n = 0;
   err_total = 0.0;
   
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
      
      DT(n,:) = [t, dt0];
      
      t = t + dt0;
      dt = min(dtmax,dt1);
      if (t+dt>tf), dt=tf-t; end
      
      err_total=err_total+dt0*MovFEM_Error_P1L2(@uexact,t,X,U,tri,tri_bf);
      fprintf('--- n = %d  t = %e dt0 = %e dt1 = %e error =  %e\n', ...
                         n,t,dt0,dt1,err_total);
      if (t>=tf-100*eps || dt < 100*eps), break; end
      
      figure(2)
      clf
      triplot(tri,X(:,1),X(:,2),'Color','r')
      axis([0 1 0 1]);
      axis square;
      drawnow;
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
      
   Ue = uexact(t,X);
   fprintf('\n N = %d  max error = %e %e\n', N, norm(Ue-U,Inf), err_total);
   
   figure(3)
   clf
   trisurf(tri,X(:,1),X(:,2),U(:,1))
   
   figure(4)
   clf
   semilogy(DT(:,1),DT(:,2));
   
   fprintf('(Nv, N) = %d %d\n', Nv, N);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uexact(t, x)

   u(:,1) = cos(2*pi*x(:,1)).*cos(3*pi*x(:,2))*exp(-0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)

   F = ut(:,1).*v(:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   G = zeros(size(x,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u;
