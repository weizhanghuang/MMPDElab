function ex1d_poisson()
%
% example for an elliptic BVP
%
% - u_xx = f in Omega = (0,1)
%
% Dirichlet BC at x = 0
% non-homogeneous Neumann BC at x = 1
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 61;

   npde = 1;
   moving_mesh = true;
   mmpde_tau = 1e-1;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   nn = 10;

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
   
% set initial conditions

   U = zeros(Nv,npde);
   
% define PDE system and BCs
   
   pdedef.bfMark = ones(Nbf,1);
   pdedef.bftype = ones(Nbf,npde);
   % neumann bc for x = 1 
   pdedef.bfMark(2,1) = 2;
   pdedef.bftype(2,npde) = zeros(1,npde); 
   
   pdedef.volumeInt = @pdedef_volumeInt;
   pdedef.boundaryInt = @pdedef_boundaryInt;
   pdedef.dirichletRes = @pdedef_dirichletRes;

% perform iteration (MP)
   
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
         Xnew = MovMesh([0,1.0],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
      else
         Xnew = X;
      end
      
      % solve physical PDEs
      
      Unew = MovFEM_bvp(U,Xnew,tri,tri_bf,pdedef,'newtons');
                        
      % update
    
      X = Xnew;
      U = Unew;
      
      figure(2)
      clf
      plot(X(:,1),U(:,1),'-o');
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);

% output
   
   err = MovFEM_Error_P1L2(@uexact,0,X,U,tri,tri_bf);
   Ue = uexact(0,X);
   fprintf('\n N = %d  max error = %e %e\n', N, norm(Ue-U,Inf), err);
   
   fprintf('(Nv, N) = %d %d\n', Nv, N);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uexact(t,x)

   u = sin(2*pi*x(:,1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)

    F = 4*pi*pi*sin(2*pi*x(:,1));
    
    F = du(:,1).*dv(:,1) - F.*v(:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   G = zeros(size(x,1),1);

   ID = find(bfMark==2);
   G(ID,1) = - 2*pi*cos(2*pi*x(ID,1)).*v(ID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u(:,1) - uexact(t,x);
