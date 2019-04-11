function ex3d_heat()
%
% example for heat equation
%
%  u_t = nabla \dot (D nabla (u)) + f
%
% on Omega = (0,1) x (0,1) x (0,1), t in (0,1], with exact solution.
%
% Dirichlet BCs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 11;

   npde = 1;
   moving_mesh = true;
   mmpde_tau = 1e-2;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   t = 0;
   tf = 1;
   dt0 = 1e-3;
   dtmax = 0.1;
   
% set the initial meshes, find the indices of the corner points and fix them

   kmax = jmax;
   lmax = jmax;
   [X,tri]=MovMesh_cube2tet(linspace(0,1,jmax),linspace(0,1,kmax),linspace(0,1,lmax));
   TR = triangulation(tri,X);
   tri_bf = freeBoundary(TR);
   Nbf = length(tri_bf);

   [Nv,d] = size(X);
   N = size(tri, 1);
   Xi_ref = X;
   
   % define nodes_fixed: corners and edges are fixed

   nodes_fixed = [];
   for i=1:Nbf
      for j=1:d
         node = tri_bf(i,j);
         if ( abs(Xi_ref(node,1)*(Xi_ref(node,1)-1))+abs(Xi_ref(node,2)*(Xi_ref(node,2)-1))<1e-8 ...
           || abs(Xi_ref(node,1)*(Xi_ref(node,1)-1))+abs(Xi_ref(node,3)*(Xi_ref(node,3)-1))<1e-8 ...
           || abs(Xi_ref(node,3)*(Xi_ref(node,3)-1))+abs(Xi_ref(node,2)*(Xi_ref(node,2)-1))<1e-8)
            nodes_fixed = [nodes_fixed, node];
         end
      end
   end
   nodes_fixed = unique(nodes_fixed);
   
% set initial conditions and compute the initial adjusted mesh

   % set the initial solution
   U = uexact(t,X);
   
   figure(1)
   clf
   plot3(X(:,1),X(:,2),X(:,3),'r+');
   view(3)
   axis tight;
   drawnow;
   
% define PDE system and BCs
   
   pdedef.bfMark = ones(Nbf,1); 
   pdedef.bftype = ones(Nbf,npde);
            
   pdedef.volumeInt = @pdedef_volumeInt;
   pdedef.boundaryInt = @pdedef_boundaryInt;
   pdedef.dirichletRes = @pdedef_dirichletRes;
   
% perform integration (MP)

   dt = dt0;
   DT = zeros(20000,2);
   err_total = 0.0;
   n = 0;
   
   tcpu = cputime;

   while true
   
      % move the mesh
      
      if (moving_mesh)
         M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
         M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
         Xnew = MovMesh([t,t+dt],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed,'ode15s');
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
      fprintf('--- n = %d  t = %e dt0 = %e dt1 = %e error = %e\n', ...
                     n,t,dt0,dt1,err_total);
      
      figure(2)
      plot3(X(:,1),X(:,2),X(:,3),'r+');
      view(3)
      axis tight;
      drawnow;
      
      if (t>=tf-100*eps || dt < 100*eps), break; end
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
      
   Ue = uexact(t,X);
   fprintf('\n N = %d  max error = %e %e\n', N, norm(Ue-U,Inf), err_total);
   
   figure(4)
   clf
   semilogy(DT(:,1),DT(:,2));
   
   fprintf('(Nv, N) = %d %d\n', Nv, N);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);
                     
%{    
   disp('ploting ...')
   figure(3)
   clf
   tetramesh(tri,X,'FaceColor','g','FaceAlpha',0.8);
   axis([0 1 0 1 0 1])
%}
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uexact(t, x)

   u(:,1) = sin(pi*x(:,1)).*sin(2*pi*x(:,2)).*sin(3*pi*x(:,3))*exp(-t);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)
    
    F = (14*pi*pi-1)*uexact(t,x);
    F=ut(:,1).*v(:)+du(:,1).*dv(:,1)+du(:,2).*dv(:,2)+du(:,3).*dv(:,3)-F.*v(:); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   G = zeros(size(x,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u(:,1)-uexact(t,x);
