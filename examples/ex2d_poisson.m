function ex2d_poisson()
%
% example for the BVP of poisson's equation
%
% - nabla( nabla(u)) = f in Omega = (0,1) x (0,1)
%
% Dirichlet BC on x = 0 & y = 0
% non-homoheneous Naumann BC on x = 1 & y = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 41;

   npde = 1;
   moving_mesh = true;
   mmpde_tau = 1e-1;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   nn = 10;
   
% set the initial meshes, find the indices of the corner points and fix them

   kmax = jmax;
   [X,tri] = MovMesh_rect2tri(linspace(0,1,jmax), linspace(0,1,kmax), 1);
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
   U = zeros(Nv,npde);
   
   figure(1)
   clf
   triplot(tri,X(:,1),X(:,2),'Color','r')
   axis([0 1 0 1]);
   axis square;

% define PDE system and BCs
   
   pdedef.bfMark = ones(Nbf,1); % for y = 0 (b1)
   Xbfm = (X(tri_bf(:,1),:)+X(tri_bf(:,2),:))*0.5;
   pdedef.bfMark(Xbfm(:,1)<1e-8) = 4; % for x = 0 (b4)
   pdedef.bfMark(Xbfm(:,1)>1-1e-8) = 2; % for x = 1 (b2)
   pdedef.bfMark(Xbfm(:,2)>1-1e-8) = 3; % for y = 1 (b3)
   
   % define boundary types
   pdedef.bftype = ones(Nbf,npde);
   % for neumann bcs:
   pdedef.bftype(pdedef.bfMark==2|pdedef.bfMark==3,npde) = 0;

   pdedef.volumeInt = @pdedef_volumeInt;
   pdedef.boundaryInt = @pdedef_boundaryInt;
   pdedef.dirichletRes = @pdedef_dirichletRes;
  
% perform integration (MP)

   tcpu = cputime;
   if (~moving_mesh)
      nn = 1;
   end

   for n=1:nn
      
      fprintf('--- n = %d\n', n);
      
      % move the mesh
      
      if (moving_mesh)

         M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
         M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
         Xnew = MovMesh([0,1],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
      else
         Xnew = X;
      end
      
      % solve physical PDEs
      
      Unew = MovFEM_bvp(U,Xnew,tri,tri_bf,pdedef,'newtons');
      %Unew = MovFEM_bvp(U,Xnew,tri,tri_bf,pdedef,'fsolve');
                        
      % update
    
      X = Xnew;
      U = Unew;
      
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
   
   figure(3)
   clf
   trisurf(tri,X(:,1),X(:,2),U(:,1))
   
   err = MovFEM_Error_P1L2(@uexact,0,X,U,tri,tri_bf);
   Ue = uexact(0, X);
   fprintf('(Nv, N, max err, L2 err) = %d %d %e %e\n',Nv,N,norm(Ue-U,Inf),err);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uexact(t, x)

   u = sin(2*pi*x(:,1)).*sin(3*pi*x(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)

    F = 13*pi*pi*sin(2*pi*x(:,1)).*sin(3*pi*x(:,2));
    F = du(:,1).*dv(:,1) + du(:,2).*dv(:,2) - F.*v(:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   G = zeros(size(x,1),1);

   ID = find(bfMark==2);
   G(ID) = -2*pi*cos(2*pi*x(ID,1)).*sin(3*pi*x(ID,2)).*v(ID);
   ID = find(bfMark==3);
   G(ID) = -3*pi*sin(2*pi*x(ID,1)).*cos(3*pi*x(ID,2)).*v(ID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = zeros(size(x,1),1);
   ID = find(bfMark==1|bfMark==4);
   Res(ID) = u(ID,1) - uexact(t,x(ID,:));
