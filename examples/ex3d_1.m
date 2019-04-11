function ex3d_1()
%
% example for adaptive mesh generation for given functions for cubic domains
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 11;

   % mmpde_int_method: 'ode15s' (default) or 'ode45'
   % mmpde_dt0: (tspan(end-1)-tspan(1))/5 or other positive numbers
   % mmpde_abstol: 1e-6 (default for ode15s), 1e-8 (default for ode15s)
   % method: 'MovMesh', 'MovMesh_XM', or 'MovMesh_X'

   mmpde_tau = 1e-2;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   mmpde_int_method = 'ode15s';
   mmpde_dt0 = [];
   mmpde_tf = 0.01;
   mmpde_abstol = [];
   method = 'MovMesh';
   nn = 5;
   
% set the initial meshes, find the indices of the corner points and fix them 

   kmax = jmax;
   lmax = jmax;
   [X,tri] = MovMesh_cube2tet(linspace(0,1,jmax),linspace(0,1,kmax),linspace(0,1,lmax));
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
   
   figure(1)
   clf
   plot3(X(:,1),X(:,2),X(:,3),'r+');
   view(3)
   axis tight;
   drawnow;

% perform integration (to generate the adaptive mesh)

   TT = zeros(nn,1);
   Ih = zeros(nn,1);
   Kmin = zeros(nn,1);

   tcpu = cputime;

   for n=1:nn
   
      % compute metric tensor
      U = funct(0, X);
      M = MovMesh_metric(U,X,tri,tri_bf,mmpde_alpha);
      M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);

      if (n==1)
         [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
      end
      
      switch (method)
      case 'MovMesh' % xi-formulation
         [Xnew,Ih(n),Kmin(n)] = MovMesh([0,mmpde_tf],Xi_ref,X,M,mmpde_tau,tri,tri_bf, ...
                        nodes_fixed,mmpde_int_method,mmpde_dt0,mmpde_abstol);
      case 'MovMesh_XM' % x-formulation
         [Xnew,Ih(n),Kmin(n)] = MovMesh_XM([0,mmpde_tf],X,M,mmpde_tau,tri,tri_bf, ...
                       nodes_fixed,mmpde_int_method,mmpde_dt0,mmpde_abstol);
      case 'MovMesh_X' % x-formulation with M = I (no adaptation)
         [Xnew,Ih(n),Kmin(n)] = MovMesh_X([0,mmpde_tf],X,mmpde_tau,tri,tri_bf, ...
                       nodes_fixed,mmpde_int_method,mmpde_dt0,mmpde_abstol);
      end
      
      fprintf('--- n = %d  Ih = %e\n', n, Ih(n));
      
      TT(n) = n;
      X = Xnew;

      figure(2)
      plot3(X(:,1),X(:,2),X(:,3),'r+');
      view(3)
      axis tight;
      drawnow;
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output

   U = funct(0,X);
   err = MovFEM_Error_P1L2(@funct,0,X,U,tri,tri_bf);
   fprintf('(Nv, N) = %d %d %e\n', Nv, N, err);
   fprintf('initial Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                          Qgeo, Qeq, Qali)
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('        Mesh quality measures (Qgeo, Qeq, Qali) = %e %e %e\n', ...
                     Qgeo, Qeq, Qali);
                     
   figure(4)
   clf
   plot(TT,Ih,'r-o');
   title('Ih')

   figure(5)
   clf
   semilogy(TT(1:end-1),abs(Ih(2:end)-Ih(1:end-1)),'-o');
   title('Ih - diff')
    
   figure(6)
   clf
   plot(TT,Kmin,'r-o');
   title('Kmin')
   
%{    
   disp('ploting ...')
   figure(3)
   clf
   tetramesh(tri,X,'FaceColor','g','FaceAlpha',0.8);
   axis([0 1 0 1 0 1])
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = funct(t, x)

u = zeros(size(x,1),1);

example = 2;

switch (example)

case 1

   u = sin(1.5*pi.*x(:,1)).*sin(1.5*pi.*x(:,2)).*sin(1.5*pi*x(:,3));

case 2

   R = 30;
   u = tanh(R*((x(:,1)-0.5).^2+(x(:,2)-0.5).^2+(x(:,3)-0.5).^2-1/16.0));
   
case 3

   R = 30;
   u = tanh(-R*(x(:,3)-0.5-0.25*sin(2*pi*x(:,1)).*sin(pi*x(:,2))));

case 4

   R = 30;
   u = tanh(-R*(x(:,3)-tanh(-R*(x(:,2)-0.5-0.25*sin(2*pi*x(:,1))))));
   
case 5

   R = 7;
   u = tanh(-R*(0.5-tanh(-R*(x(:,2)-0.5-0.25*sin(2*pi* x(:,1))))));
end


