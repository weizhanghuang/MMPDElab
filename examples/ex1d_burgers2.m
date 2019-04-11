function ex1d_burgers2()
%
% example for burgers equation
%
% u_t = epsilon u_xx - u u_x
%
% on Omega = (0,1), t in (0,1]
%
% Dirichlet BC on x = 0 & x = 1
%
% smooth initial solution, exact solution not available
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global epsilon;

if (~isdeployed)
  addpath('../src_MMPDElab');
end

movie_draw = false;
if (movie_draw)
   iiout = 0;
   nout = 1;
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
   dtmax = 0.01;
   
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
   U = uinitial(X);
   
   % generate initial adjusted mesh
   if (moving_mesh)
       for n=1:5  
          M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
          M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
          X = MovMesh([0,1],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
          U = uinitial(X);
             
          figure(1)
          plot(X,U,'-o');
          hold on
          plot(X,ones(Nv,1)*(-0.1),'+r');
          hold off
          axis([0 1 -1 1.5]);
          drawnow;
       end
   end

% define PDE system and BCs

   % all bcs are dirichlet so no need for marking boundary segments
   pdedef.bfMark = ones(Nbf,1);
   pdedef.bftype = ones(Nbf,npde);

   pdedef.volumeInt = @pdedef_volumeInt;
   pdedef.boundaryInt = @pdedef_boundaryInt;
   pdedef.dirichletRes = @pdedef_dirichletRes;
   
% perform integration (MP)

   DT = zeros(2000,2);
   tcpu = cputime;
   n = 0;
   dt = dt0;

   while true
      
      % move the mesh
      
      if (moving_mesh)
         M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
         M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
         M = Matrix_ceil(M,1e4);
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
      
      figure(2)
      plot(X,U,'-o');
      hold on
      plot(X,ones(Nv,1)*(-0.8),'+r');
      hold off
      axis([0 1 -1 1.5]);
      drawnow;
      
if (movie_draw & mod(n,nout)==0)
      iiout = iiout + 1;
      figure(18)
      clf
      plot(X,U,'r-o');
      hold on
      plot(X,-0.8*ones(Nv,1),'r+');
      hold off
      axis([0 1 -1 1.5]);
      title(['t = ',num2str(t)])
      drawnow;
      pngName1 = strcat('mesh',num2str(iiout),'.png');
      print(18,'-dpng',pngName1);
end
      
      fprintf('--- n = %d  t = %e dt0 = %e dt1 = %e\n',n,t,dt0,dt1);
      
      if (t>=tf-100*eps || dt < 100*eps), break; end
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);

% output
         
   figure(4)
   clf
   semilogy(DT(:,1),DT(:,2));
   title('dt')
   
   fprintf('(Nv, N, max err, L2 err) = %d %d\n',Nv,N);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uinitial(x)

   u(:,1) = sin(2*pi*x(:,1)) + 0.5*sin(pi*x(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)

   global epsilon
    
   F = ut(:,1).*v(:) + epsilon*du(:,1).*dv(:,1) + u(:,1).*du(:,1).*v(:); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   G = zeros(size(x,1),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u(:,1)-uinitial(x);


