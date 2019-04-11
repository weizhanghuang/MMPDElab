function ex2d_combustion()
%
% example for a combusition model consisting of two PDEs
% on a complex domain and subject to Dirichlet and Neumann BCs.
% exact solution is not available.
%
% the description of the model can be found in: Jochen Froehlich and Jens Lang,
% Two-dimensional cascadic finite element computations of combustion problems,
% Comput. Methods Appl. Mech. Engrg. 158 (1998) 255-267.
%
% this example is provided by Jens Lang and Alf Gerisch (TU Darmstadt, Germany)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   kmax = 3;

   npde = 2;
   moving_mesh = true;
   mmpde_tau = 1e-3;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   t = 0;
   tf = 60;
   dt0 = 1e-4;
   dtmax = 0.2;

% set the initial meshes, find the indices of the corner points and fix them

   jmax = 4*kmax;
   [X1,tri1] = MovMesh_rect2tri(linspace(0,15,jmax+1), ...
                                linspace(0,16,4*kmax+1),1);
   [X2,tri2] = MovMesh_rect2tri(linspace(15,30,jmax+1), ...
                                linspace(4,12,2*kmax+1),1);
   [X,tri] = MovMesh_MeshMerge(X1,tri1,X2,tri2);
   [X3,tri3] = MovMesh_rect2tri(linspace(30,60,2*jmax+1), ...
                                linspace(0,16,4*kmax+1),1);
   [X,tri] = MovMesh_MeshMerge(X,tri,X3,tri3);
   TR = triangulation(tri,X);
   tri_bf = freeBoundary(TR);
   Nbf = length(tri_bf);

   [Nv,d] = size(X);
   N = size(tri, 1);
   Xi_ref = X;

   % find the indices of the corner points and fix them
   corners = [0,0;15,0;15,4;30,4;30,0;60,0;60,16;30,16;30,12;15,12;15,16;0,16];
   [~,nodes_fixed] = ismembertol(corners,Xi_ref,1e-10,'ByRows',true);
   % nodes_fixed = unique(tri_bf); % for fixing all boundary nodes
   
% set initial conditions

   % set the initial solution
   U = uinitial(X);
   
   % generate initial adjusted mesh
   if (moving_mesh)   
      for n=1:5 
         M = MovMesh_metric(U(:,1),X,tri,tri_bf,mmpde_alpha);
         M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
         Xnew = MovMesh([0,1e-2],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
         X = Xnew;
         U = uinitial(X); 
      
         figure(1)
         clf
         triplot(tri,X(:,1),X(:,2),'Color','r')
         axis([0 60 0 16]);
         axis equal; 
         drawnow;
      end
   end
   
% define PDE system and BCs

   pdedef.bfMark = ones(Nbf,1);  
   Xbfm = (X(tri_bf(:,1),:)+X(tri_bf(:,2),:))*0.5;
   pdedef.bfMark(Xbfm(:,1) < 1e-8) = 2;
   pdedef.bfMark(abs(Xbfm(:,1)-15) < 1e-8) = 3;
   pdedef.bfMark(abs(Xbfm(:,1)-30) < 1e-8) = 3;
   pdedef.bfMark((abs(Xbfm(:,2)-4) < 1e-8) & ...
       (Xbfm(:,1) > 15 & Xbfm(:,1) < 30)) = 3; 
   pdedef.bfMark((abs(Xbfm(:,2)-12) < 1e-8) & ...
       (Xbfm(:,1) > 15 & Xbfm(:,1) < 30)) = 3; 

   pdedef.bftype = ones(Nbf,npde);
   pdedef.bftype(pdedef.bfMark==1|pdedef.bfMark==3,:) = 0; 
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
         M = MovMesh_metric(U(:,1),X,tri,tri_bf,mmpde_alpha);
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
      
      fprintf('--- n = %d  t = %e dt0 = %e dt1 = %e\n',n,t,dt0,dt1);
      
      figure(2)
      clf
      triplot(tri,X(:,1),X(:,2),'Color','r')
      axis([0 60 0 16]);
      axis equal; 
      drawnow;
      
      figure(3)
      clf
      subplot(2,1,1);
      trisurf(tri,X(:,1),X(:,2),U(:,1));
      title(['temperature at t = ' num2str(t)]);
      axis([0 60 0 16]);
      axis equal;
      subplot(2,1,2);
      trisurf(tri,X(:,1),X(:,2),U(:,2));
      title(['fuel at t = ' num2str(t)]);
      axis([0 60 0 16]);
      axis equal;
      drawnow;
      
      if (t>=tf-100*eps || dt < 100*eps), break; end
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
   
   figure(5)
   clf
   semilogy(DT(:,1),DT(:,2));
   
   fprintf('(Nv, N) = %d %d\n', Nv, N);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);
                     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uinitial(x)

   x0 = 7.5;
   Le = 1;
   u = zeros(size(x,1),2);
   
   u(:,1) = 1*(x(:,1)<=x0) + ( exp(x0-x(:,1)) ).*(x(:,1)>x0);
   u(:,2) = 0*(x(:,1)<=x0) + ( 1-exp(Le*(x0-x(:,1))) ).*(x(:,1)>x0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)
    
    beta = 10;
    alpha = 0.8;
    Le = 1;
    
    w = beta^2/(2*Le)*u(:,2).*exp(-beta*(1-u(:,1))./(1-alpha*(1-u(:,1))));
    
    if (ipde==1)
       F = ut(:,1).*v+du(:,1).*dv(:,1)+du(:,2).*dv(:,2) - w.*v;
    else
       F = ut(:,2).*v+(du(:,3).*dv(:,1)+du(:,4).*dv(:,2))/Le + w.*v;
    end 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   k = 0.1;
   G = zeros(size(x,1),1);
   
   if ipde==1
      ID = find(bfMark==3);
      G(ID) = k*u(ID,1).*v(ID);
   end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = zeros(size(x,1),1);
   ID = find(bfMark==2);
   if (ipde==1)
      Res(ID) = u(ID,1)-1;
   else
      Res(ID) = u(ID,2)-0;
   end
