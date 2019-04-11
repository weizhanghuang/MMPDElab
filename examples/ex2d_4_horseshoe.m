function ex2d_4_horseshoe()
%
% example for horseshoe shape -- smooth the algebraic mesh
%
% in this example, the computational and physical meshes are different.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (~isdeployed)
  addpath('../src_MMPDElab');
end

% set the basic parameters

   jmax = 31;

   mmpde_tau = 1e-2;
   mmpde_ncycles = 3;
   mmpde_tf = 0.1;
   nn = 10;
   
% set the initial meshes, find the indices of the corner points and fix them

   kmax = jmax;
   [X,tri] = MovMesh_rect2tri(linspace(0,1,jmax),linspace(0,1,kmax),1);
   TR = triangulation(tri,X);
   tri_bf = freeBoundary(TR);
   Nbf = length(tri_bf);
   
   [Nv,d] = size(X);
   N = size(tri, 1);
   Xi_ref = X;
   
   R = 4.5;
   X(:,1) = - (1+Xi_ref(:,2)).*cos(pi*Xi_ref(:,1));
   X(:,2) = (1+(2*R-1)*Xi_ref(:,2)).*sin(pi*Xi_ref(:,1));
   
   % find the indices of the corner points and fix them
   corners = [0, 0; 1, 0; 1, 1; 0, 1]; % for xi-mesh
   [~,nodes_fixed] = ismembertol(corners,Xi_ref,1e-10,'ByRows',true);
   % fix nodes on the top boundary are fixed
   ID = find(Xi_ref(:,2)>1-1e-8);
   nodes_fixed = [nodes_fixed;ID];
   nodes_fixed = unique(nodes_fixed);

   figure(1)
   clf
   triplot(tri,X(:,1),X(:,2),'Color','r')
   axis([-2, 2, 0, 2*R])
   axis equal;
   drawnow;
      
% perform integration (to generate the adaptive mesh)

   TT = zeros(nn,1);
   Ih = zeros(nn,1);
   Kmin = zeros(nn,1);
   M = repmat(reshape(eye(d,d),1,[]),Nv,1);

   tcpu = cputime;

   for n=1:nn
   
      if (n==1)
         [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
      end

      % compute metrix tensor
%{
      temp = 1.0 + 1.0./(exp(4*sqrt(X(:,1).^2 + (X(:,2)-2*R).^2))-1+1/sqrt(N));
      M = repmat(reshape(eye(d,d),1,[]),Nv,1);
      M = bsxfun(@times,M,temp); 
      M = MovMesh_metric_smoothing(M,mmpde_ncycles,X,tri);
%}
      
      [Xnew,Ih(n),Kmin(n)] = MovMesh([0,mmpde_tf],Xi_ref,X,M,mmpde_tau,tri,tri_bf,nodes_fixed);
      
      fprintf('--- n = %d  Ih = %e\n', n, Ih(n));
            
      TT(n) = n;
      X = Xnew;
       
      figure(2)
      triplot(tri,X(:,1),X(:,2),'Color','r')
      axis([-2, 2, 0, 2*R])
      axis equal;
      drawnow;
      
   end
   
   tcpu = cputime-tcpu;
   fprintf('\n --- total elapsed cpu time = %e \n\n', tcpu);

% output
   
   fprintf('(Nv, N, error) = %d %d\n', Nv, N);
   fprintf('initial mesh (Qgeo, Qeq, Qali) = %e %e %e\n',Qgeo,Qeq,Qali)
   [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M);
   fprintf('final   Mesh (Qgeo, Qeq, Qali) = %e %e %e\n',Qgeo,Qeq,Qali)
   
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
