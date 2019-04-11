function ex2d_1()
%
% example for adaptive mesh generation for given functions on square domains
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 31;

   % mmpde_int_method: 'ode15s' (default) or 'ode45'
   % mmpde_dt0: (tspan(end-1)-tspan(1))/5 or other positive numbers
   % mmpde_abstol: 1e-6 (default for ode15s), 1e-8 (default for ode15s)
   % method: 'MovMesh', 'MovMesh_XM', or 'MovMesh_X'

   mmpde_tau = 1e-2;
   mmpde_ncycles = 3;
   mmpde_alpha = [];
   mmpde_int_method = [];
   mmpde_dt0 = [];
   mmpde_tf = 0.1;
   mmpde_abstol = [];
   method = 'MovMesh';
   nn = 10;
   
% set the initial meshes, find the indices of the corner points and fix them

   kmax = jmax;
   xmin = 0;
   xmax = 1;
   ymin = 0;
   ymax = 1;
   [X,tri] = MovMesh_rect2tri(linspace(0,1,jmax),linspace(0,1,kmax),1);
   TR = triangulation(tri,X);
   tri_bf = freeBoundary(TR);
   Nbf = length(tri_bf);
   
   [Nv,d] = size(X);
   N = size(tri, 1);
   Xi_ref = X;
   
   % perturb the mesh
%{      
   h = (1/N)^(1/d);
   dX = (2*rand(Nv,d)-1)*0.5*h*0.8;
   dX(tri_bf,:) = 0;
   X = X + dX;
%}

   % find the indices of the corner points and fix them
   corners = [xmin, ymin; xmax, ymin; xmax, ymax; xmin, ymax];
   [~,nodes_fixed] = ismembertol(corners,Xi_ref,1e-10,'ByRows',true);
   % nodes_fixed = unique(tri_bf); % for fixing all boundary nodes

   figure(1)
   clf
   triplot(tri,X(:,1),X(:,2),'Color','r')
   axis square;
   drawnow;
      
% perform integration (to generate the adaptive mesh)
   
   TT = zeros(nn,1);
   Ih = zeros(nn,1);
   Kmin = zeros(nn,1);

   tcpu = cputime;

   for n=1:nn

      % compute metric tensor
      U = uexact(0, X);
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
      triplot(tri,X(:,1),X(:,2),'Color','r')
      axis square;
      drawnow;
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);


% output
   
   U = uexact(0,X);
   err = MovFEM_Error_P1L2(@uexact,0,X,U,tri,tri_bf);
   fprintf('(Nv, N, error) = %d %d %e\n', Nv, N, err);
   fprintf('initial mesh (Qgeo, Qeq, Qali) = %e %e %e\n',Qgeo,Qeq,Qali)
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('final   Mesh (Qgeo, Qeq, Qali) = %e %e %e\n',Qgeo,Qeq,Qali)

   figure(3)
   clf
   trisurf(tri,X(:,1),X(:,2),U(:,1),'EdgeColor','none','LineStyle','none', ...
          'FaceLighting','phong')
   view(2)
   axis square;
   colorbar
   
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
   semilogy(TT,Kmin,'r-o');
   title('Kmin')
    
   figure(7)
   clf
   trisurf(tri,X(:,1),X(:,2),U(:,1),'EdgeColor','none','LineStyle', ...
           'none','FaceColor','interp');
   title('U(:,1)')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u = uexact(t, x)

u = zeros(size(x,1),1);

example = 3;

switch (example)

case 1

   u = sin(2*pi.*x(:,1)).*sin(3*pi.*x(:,2));

case 2

   R = 30.0;
   u = tanh(R*((x(:,1)-0.5).^2+(x(:,2)-0.5).^2-1.0/16.0));
   
case 3

   R = 30.0;
   u = tanh( -R * ( x(:,2) - 0.5 - 0.25*sin(2*pi*x(:,1)) ) );
  
case 4

   R = 30;
   u = tanh(R*2*x(:,2))-tanh(0.5*R*4*x(:,1)-R*2*x(:,2)-R);
   
case 5

   R = 30.0;
   xx = -2 + 4*x(:,1);
   yy = -2 + 4*x(:,2);
   u = tanh(R*(xx.^2+yy.^2-1.0/8)) ...
     + tanh(R*((xx-0.5).^2+(yy-0.5).^2-1.0/8)) ...
     + tanh(R*((xx-0.5).^2+(yy+0.5).^2-1.0/8)) ...
     + tanh(R*((xx+0.5).^2+(yy-0.5).^2-1.0/8)) ...
     + tanh(R*((xx+0.5).^2+(yy+0.5).^2-1.0/8));
end


