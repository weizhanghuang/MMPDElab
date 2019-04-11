function ex2d_poisson_Lshape()
%
% example for an elliptic BVP
%
%    - beta^2 * nabla( nabla(u)) + u = 0
%
% in L-shaped domain, subject to
%
% Dirichlet BC  u = 1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

global beta

% set the basic parameters

   jmax = 21;

   beta = 2^(-5);
   npde = 1;
   moving_mesh = true;
   mmpde_tau = 1e-2;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   nn = 10;
   
% set the initial meshes, find the indices of the corner points and fix them

   kmax = jmax;
   [X1,tri1] = MovMesh_rect2tri(linspace(0,1,jmax),linspace(0,1,kmax),1);
   [X2,tri2] = MovMesh_rect2tri(linspace(1,2,jmax),linspace(0,1,kmax),1);
   [X3,tri3] = MovMesh_rect2tri(linspace(0,1,jmax),linspace(1,1.5,(kmax-1)/2+1),1);
   [X,tri] = MovMesh_MeshMerge(X1,tri1,X2,tri2);
   [X,tri] = MovMesh_MeshMerge(X,tri,X3,tri3);
   TR = triangulation(tri,X);
   tri_bf = freeBoundary(TR);
   Nbf = length(tri_bf);

   [Nv,d] = size(X);
   N = size(tri, 1);
   Xi_ref = X;

   % find the indices of the corner points and fix them
   corners = [0,0;2,0;2,1;1,1;1,1.5;0,1.5];
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
   
   pdedef.bfMark = ones(Nbf,1);
   pdedef.bftype = ones(Nbf,npde);
   
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
         M_max = max(max(abs(M)));
         temp = 1.0 + 1.0./(exp(4*sqrt((X(:,1)-1).^2 + (X(:,2)-1).^2))-1+1/M_max);
         M = M + bsxfun(@times,repmat(reshape(eye(d,d),1,[]),Nv,1),temp); 
         M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
         Xnew = MovMesh([0,1],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
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
      triplot(tri,X(:,1),X(:,2),'Color','r')
      axis([0 2 0 1.5]);
      axis square;
      drawnow;
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
      
   figure(3)
   clf
   trisurf(tri,X(:,1),X(:,2),U(:,1))
   
   figure(4)
   clf
   trisurf(tri,X(:,1),X(:,2),-beta*log(abs(U(:,1))))
   
   fprintf('(Nv, N) = %d %d\n', Nv, N);
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function F = pdedef_volumeInt(du, u, ut, dv, v, x, t, ipde)

global beta;

    d11 = 1*beta^2;
    d12 = 0*beta^2;
    d22 = 1*beta^2;

    F = d11*du(:,1).*dv(:,1) + d12*du(:,2).*dv(:,1) + d12*du(:,1).*dv(:,2) ...
      + d22*du(:,2).*dv(:,2) + u.*v; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function G = pdedef_boundaryInt(du, u, v, x, t, ipde, bfMark)

   G = zeros(size(x,1),1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Res = pdedef_dirichletRes(u, x, t, ipde, bfMark)

   Res = u(:,1) - 1;

