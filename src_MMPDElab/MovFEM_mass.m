function B = MovFEM_mass(t,u,udot,npde,X,Xdot,tri,tri_bf,tn,pdedef)
%
% usage: B = MovFEM_mass(t,u,udot,npde,X,Xdot,tri,tri_bf,tn,pdedef)
%
% this function computes F_udot(t,u,udot) (or dfdyt) for the system
% F(t,u,udot) = 0 resulting from the P1 conforming fem discretization
% of IBVPs described in MovFEM(). it uses finite difference approximation.
%
% t:    current time.
% u:    u = (u_1, u_2, ..., u_npde) the unknown variables, of size Nv*npde-by-1,
%       u(npde*(j-1)+k), j=1:Nv, k=1:npde.
% udot: (material) time derivatives of the unknown variables,
%       of size Nv*npde-by-1, udot(npde*(j-1)+k), j=1:Nv, k=1:npde.
% npde: number of the components of the physical solution.
% X:    coordinates of mesh vertices at time t_n, of size Nv-by-d.
% Xdot: mesh velocity, of size Nv-by-d.
% tri:  the connectivity for all meshes, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets for all meshes, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% tn:   beginning time at the current step.
% pdedef: a structure used to define the PDE system in weak form; see MovFEM().
% B:    (output) mass matrix (sparse) of dimension Nv*npde-by-Nv*npde.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

% compute basic parameters

    [Nv,d] = size(X);
    N = size(tri,1);
    Nbf = size(tri_bf,1);

    U = reshape(u,npde,Nv)';
    Udot = reshape(udot,npde,Nv)';
    x = X + (t-tn)*Xdot;

% integration nodes and weights

   [lambda,nn] = MovFEM_int_weights(d);
   
% compute edge matrix and Phi (gradients of basis functions).
% Phi(:,(j-1)*d+(1:d)) = \nabla phi_j

    [~, Einv, detE, ~] = Matrix_edge(x, tri);
    Phi = zeros(N,d*(d+1));
    for j=1:d
    for k=1:d
        Phi(:,j) = Phi(:,j) - Einv(:,(j-1)*d+k);
    end
    end
    for j=1:d
    for k=1:d
        Phi(:,d+(j-1)*d+k) = Einv(:,(k-1)*d+j);
    end
    end
    detE = abs(detE)/factorial(d);
   
% compute dU on elements: dU(:,(i-1)*d+(1:d)) = \nabla u_i
   
    dU = zeros(N,d*npde);
    for i=1:npde
    for j=1:(d+1)
    for k=1:d
        dU(:,d*(i-1)+k) = dU(:,d*(i-1)+k) ...
                        + Phi(:,(j-1)*d+k).*U(tri(:,j),i);
    end
    end
    end
   
% compute rhs on elements (volume integrals)

    FF = zeros(N,npde*(d+1));
    for n = 1:nn
        xK = zeros(N,d);
        xdotK = zeros(N,d);
        uK = zeros(N,npde);
        udotK = zeros(N,npde);
        for j=1:(d+1)
            xK = xK + x(tri(:,j),:)*lambda(j,n);
            xdotK = xdotK + Xdot(tri(:,j),:)*lambda(j,n);
            uK = uK + U(tri(:,j),:)*lambda(j,n);
            udotK = udotK + Udot(tri(:,j),:)*lambda(j,n);
        end
        for i=1:npde
        for k=1:d
            udotK(:,i) = udotK(:,i)-dU(:,d*(i-1)+k).*xdotK(:,k);
        end
        end
        for j=1:(d+1)
            phiK = repmat(lambda(j,n),N,1);
            dphiK = Phi(:,(d*(j-1)+1):(d*j));
            for i=1:npde
               FF(:,npde*(j-1)+i) = FF(:,npde*(j-1)+i) ...
             + pdedef.volumeInt(dU,uK,udotK,dphiK,phiK,xK,t,i)*lambda(d+2,n);        
            end
        end
    end
    FF = bsxfun(@times,FF,detE);

% compute jacobian with respect to udot term

len_b = 0;
for ipde=1:npde
   nodes = unique(tri_bf(pdedef.bftype(:,ipde)==1,:));
   len_b = max(len_b,length(nodes));
end

II = zeros(N*npde*npde*(d+1)*(d+1)+Nbf*npde*npde*d*d+len_b*len_b*npde*npde,1);
JJ = zeros(N*npde*npde*(d+1)*(d+1)+Nbf*npde*npde*d*d+len_b*len_b*npde*npde,1);
B_nz = zeros(N*npde*npde*(d+1)*(d+1)+Nbf*npde*npde*d*d+len_b*len_b*npde*npde,1);

ic = 0;
Indx = zeros(N,npde*(d+1));

UUdot = zeros(N,npde,d+1);
for j=1:(d+1)
   UUdot(:,:,j)= Udot(tri(:,j),:);
end

for j1=1:(d+1)
for k1=1:npde
    
    h = (UUdot(:,k1,j1)+max(abs(UUdot(:,k1,j1)),1)*sqrt(eps))-UUdot(:,k1,j1);
    ss = (sign(UUdot(:,k1,j1))>= 0);
    h = (ss-(~ss)).*abs(h);
    
    UUdot(:,k1,j1) = UUdot(:,k1,j1) + h; 

    FF1 = zeros(N,npde*(d+1));
    for n = 1:nn
        xK = zeros(N,d);
        xdotK = zeros(N,d);
        uK = zeros(N,npde);
        udotK = zeros(N,npde);
        for j=1:(d+1)
            xK = xK + x(tri(:,j),:)*lambda(j,n);
            xdotK = xdotK + Xdot(tri(:,j),:)*lambda(j,n);
            uK = uK + U(tri(:,j),:)*lambda(j,n);
            udotK = udotK + UUdot(:,:,j)*lambda(j,n);
        end
        for i=1:npde
        for k=1:d
            udotK(:,i) = udotK(:,i)-dU(:,d*(i-1)+k).*xdotK(:,k);
        end
        end
        for j=1:(d+1)
            phiK = repmat(lambda(j,n),N,1);
            dphiK = Phi(:,(d*(j-1)+1):(d*j));
            for i=1:npde
               FF1(:,npde*(j-1)+i) = FF1(:,npde*(j-1)+i) ...
              + pdedef.volumeInt(dU,uK,udotK,dphiK,phiK,xK,t,i)*lambda(d+2,n);        
            end
        end
    end
    FF1 = bsxfun(@times,FF1,detE);  
    UUdot(:,k1,j1) = UUdot(:,k1,j1) - h; 
   
    h = 1.0./h;
    FF1 = FF1-FF;
    FF1 = bsxfun(@times,FF1,h);
    B_nz(ic+1:ic+N*npde*(d+1)) = reshape(FF1, N*npde*(d+1), 1);
    for j=1:(npde*(d+1))
        Indx(:,j) = npde*(tri(:,j1)-1)+k1;
    end
    JJ(ic+1:ic+N*npde*(d+1)) = reshape(Indx, N*npde*(d+1), 1);
    for j=1:(d+1)
    for k=1:npde
        Indx(:,(j-1)*npde+k) = npde*(tri(:,j)-1)+k;
    end
    end
    II(ic+1:ic+N*npde*(d+1)) = reshape(Indx, N*npde*(d+1), 1);
    ic = ic + N*npde*(d+1);
end
end

% save to sparse matrix

    B = sparse(II(1:ic), JJ(1:ic), B_nz(1:ic), Nv*npde, Nv*npde);
    
% for dirichlet boundary conditions
    % for rows
    B = B';
    for ipde=1:npde
        for j=1:d
            nodes = tri_bf(pdedef.bftype(:,ipde)==1,j);
            B(:,npde*(nodes-1)+ipde) = 0;
        end
    end
    B = B';
    
% for lumping

   lumping = false;
   if (lumping)
      B = B';
      bb = sum(B)';
      B = spdiags(bb,0,Nv*npde,Nv*npde);
   end
    
% end of MovFEM_mass()
