function ex2d_burgers_2level()
%
% example for burgers equation (with 2-level moving mesh strategy)
%
%        u_t = epsilon \Delta u - u u_x - u u_v
%
% in Omega = (-0.5,1) x (-0.5,1), t in (0,2]
%
% u = 0 on \partial \Omega
%
% Dirichlet BC
% exact solution is not available
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global epsilon;

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 31;
   Level = 1;

   epsilon = 5e-3;
   npde = 1;
   moving_mesh = true;
   mmpde_tau = 1e-3;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   t = 0;
   tf = 2;
   dt0 = 1e-2;
   dtmax = 0.025;

% set the initial meshes, find the indices of the corner points and fix them

   kmax = jmax;
   a = -0.5;
   b = 1.0;
   [X,tri] = MovMesh_rect2tri(linspace(a,b,jmax), linspace(a,b,kmax), 1);
   TR = triangulation(tri,X);
   tri_bf = freeBoundary(TR);
   Nbf = length(tri_bf);

   [Nv,d] = size(X);
   N = size(tri, 1);
   Xi_ref = X;

   % find the indices of the corner points and fix them
   corners = [a, a; b, a; b, b; a, b];
   [~,nodes_fixed] = ismembertol(corners,Xi_ref,1e-10,'ByRows',true);
   % nodes_fixed = unique(tri_bf); % for fixing all boundary nodes
   
% compute the fine mesh and its relation to the coarse mesh

   XC = X;
   TriC = tri;
   TriC_bf = tri_bf;
   [X,tri,tri_parent] = MovMesh_MeshUniformRefine(XC,TriC,Level);
   TR = triangulation(tri,X);
   tri_bf = freeBoundary(TR);
   NvC = Nv;
   NC = N;
   Nv = size(X,1);
   N = size(tri,1);
   
% set initial conditions and compute the initial adjusted mesh

   % set the initial solution
   U = uinitial(t,X);
   
   if (moving_mesh)   
      for n=1:5    
         M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
         % project the metric tensor from fine mesh to coarse mesh
         MC = MovMesh_metric_F2C(M,tri,tri_parent,TriC);
         MC = MovMesh_metric_smoothing(MC,mmpde_ncycles,XC,TriC);
         % move the coarse mesh
         XC = MovMesh([0,1],Xi_ref,XC,MC,mmpde_tau,TriC,TriC_bf,nodes_fixed);
         % compute the new fine mesh
         [X,~,~] = MovMesh_MeshUniformRefine(XC,TriC,Level);
         U = uinitial(t,X);
      
         figure(1)
         triplot(tri,X(:,1),X(:,2),'Color','r')
         axis([a b a b]);
         axis square;
         drawnow;
      end
   end
   
% define PDE system and BCs
   
   % no need to mark the boundary segments since dirichlet bcs for all
   % boundary segments.
   
   pdedef.bfMark = ones(Nbf,1);    
   pdedef.bftype = ones(Nbf,npde);
   
   pdedef.volumeInt = @pdedef_volumeInt;
   pdedef.boundaryInt = @pdedef_boundaryInt;
   pdedef.dirichletRes = @pdedef_dirichletRes;

% perform integration (MP)

   dt = dt0;
   DT = zeros(20000,2);   
   tcpu = cputime;
   n = 0;

   while true
      
      % move the mesh
      
      if (moving_mesh)
         M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
         MC = MovMesh_metric_F2C(M,tri,tri_parent,TriC);
         MC = MovMesh_metric_smoothing(MC,mmpde_ncycles,XC,TriC);
         XC = MovMesh([t,t+dt],Xi_ref,XC,MC,mmpde_tau,TriC,TriC_bf,nodes_fixed);
         [Xnew,~,~] = MovMesh_MeshUniformRefine(XC,TriC,Level);
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
      
      figure(2)
      triplot(tri,X(:,1),X(:,2),'Color','r')
      axis([a b a b]);
      axis square;
      drawnow;
      
      fprintf('--- n = %d  t = %e dt0 = %e dt1 = %e\n',n,t,dt0,dt1);
      
      if (t>=tf-100*eps || dt < 100*eps), break; end
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
   
   figure(3)
   clf
   trisurf(tri,X(:,1),X(:,2),U(:,1));
   axis([a b a b]);
   axis square;
   colorbar
   
   figure(4)
   clf
   semilogy(DT(:,1),DT(:,2));
   
   figure(5)
   clf
   [xx,yy] = meshgrid(linspace(a,b,10*(jmax-1)+1),linspace(a,b,10*(kmax-1)+1));
   uu = griddata(X(:,1),X(:,2),U(:,1),xx,yy);
   contourf(xx,yy,uu,10);
   axis square;
   colorbar
   
   fprintf('(Nv, N, NvC, NC) = %d %d %d %d\n', Nv, N, NvC, NC);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u = uinitial(t, x)

global epsilon;

   c = - log(10^(-16));
   u = exp(-c*(x(:,1).^2+x(:,2).^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)

    global epsilon;
   
    F = ut(:,1).*v(:) + epsilon*(du(:,1).*dv(:,1) + du(:,2).*dv(:,2)) ...
        + u(:).*(du(:,1).*v(:) + du(:,2).*v(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   G = zeros(size(x,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u(:,1) - 0.0;


