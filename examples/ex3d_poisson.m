function ex3d_poisson()
%
% example for an elliptic BVP
%
% - (u_xx + u_yy + u_zz) = f
%
% in Omega = (0,1) x (0,1) x (0,1). 
% nonhomogeneous Neumann BC on X = 1 and Dirichlet BCs for other part of boundary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 11;

   npde = 1;
   moving_mesh = true;
   mmpde_tau = 1e-1;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   nn = 5;

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

   % define nodes_fixed: corners and boundary edges are fixed
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

% set initial conditions

   % set the initial solution
   U = zeros(Nv,npde);
   
   figure(1)
   clf
   plot3(X(:,1),X(:,2),X(:,3),'r+');
   view(3)
   axis tight;
   drawnow;

% define PDE system and BCs

   pdedef.bfMark = ones(Nbf,1);
   Xbfm = (X(tri_bf(:,1),:)+X(tri_bf(:,2),:)+X(tri_bf(:,3),:))/3;
   pdedef.bfMark(Xbfm(:,1)>1-1e-8) = 2; % for x=1
   pdedef.bftype = ones(Nbf,npde);
   pdedef.bftype(pdedef.bfMark==2,npde) = 0; % neumann bc for x=1
   
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
         Xnew = MovMesh([0,1],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed,'ode15s');
      else
         Xnew = X;
      end
      
      % solve physical PDEs
      
      Unew = MovFEM_bvp(U,Xnew,tri,tri_bf,pdedef,'newtons');
                        
      % update
    
      X = Xnew;
      U = Unew;
      
      figure(2)
      plot3(X(:,1),X(:,2),X(:,3),'r+');
      view(3)
      axis tight;
      drawnow;
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
      
   err = MovFEM_Error_P1L2(@uexact,0,X,U,tri,tri_bf);
   Ue = uexact(0, X);
   fprintf('\n N = %d  max error = %e %e\n', N, norm(Ue-U,Inf), err);
   
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

    u = sin(2*pi*x(:,1)).*sin(3*pi*x(:,2)).*sin(pi*x(:,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)
    
    F = 14*pi^2*sin(2*pi*x(:,1)).*sin(3*pi*x(:,2)).*sin(pi*x(:,3));
    F = du(:,1).*dv(:,1)+du(:,2).*dv(:,2)+du(:,3).*dv(:,3)-F.*v(:); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   G = zeros(size(x,1),1);
   ID = find(bfMark==2);
   G(ID) = -2*pi*sin(3*pi*x(ID,2)).*sin(pi*x(ID,3)).*v(ID);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u(:,1) - uexact(t,x);

