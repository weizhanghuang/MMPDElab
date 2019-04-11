function JAC = MovFEM_jac(t,u,udot,npde,Xn,Xdot,tri,tri_bf,tn,pdedef)
%
% usage: JAC = MovFEM_jac(t,u,udot,npde,Xn,Xdot,tri,tri_bf,tn,pdedef)
%
% this function computes F_u(t,u,udot) for the system F(t,u,udot) = 0
% resulting from the P1 conforming fem discretization of IBVPs described
% in MovFEM(). it uses finite difference approximation.
%
% t:    current time.
% u:    u = (u_1, u_2, ..., u_npde) the unknown variables, of size Nv*npde-by-1,
%       u(npde*(j-1)+k), j=1:Nv, k=1:npde.
% udot: (material) time derivatives of the unknown variables,
%       of size Nv*npde-by-1, udot(npde*(j-1)+k), j=1:Nv, k=1:npde.
% npde: number of the components of the physical solution.
% Xn:   coordinates of mesh vertices at time t_n, of size Nv-by-d.
% Xdot: mesh velocity, of size Nv-by-d.
% tri:  the connectivity for all meshes, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets for all meshes, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% tn:   beginning time at the current step.
% pdedef: a structure used to define the PDE system in the weak form;
%       see MovFEM().
% JAC:  (output) jacobian matrix (sparse) of dimension Nv*npde-by-Nv*npde.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

% compute basic parameters

    [Nv,d] = size(Xn);
    N = size(tri,1);
    Nbf = size(tri_bf,1);

    U = reshape(u,npde,Nv)';
    Udot = reshape(udot,npde,Nv)';
    X = Xn + (t-tn)*Xdot;

% integration nodes and weights

   [lambda,nn] = MovFEM_int_weights(d);
   
% compute edge matrix and Phi (gradients of basis functions).
% Phi(:,(j-1)*d+(1:d)) = \nabla phi_j

    [~, Einv, detE, ~] = Matrix_edge(X, tri);
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
                        + Phi(:,(j-1)*d+k).*U((tri(:,j)),i);
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
            xK = xK + X(tri(:,j),:)*lambda(j,n);
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

% compute jacobian

len_b = 0;
for ipde=1:npde
   nodes = unique(tri_bf(pdedef.bftype(:,ipde)==1,:));
   len_b = max(len_b,length(nodes));
end

II = zeros(N*npde*npde*(d+1)*(d+1)+Nbf*npde*npde*d*d+len_b*len_b*npde*npde,1);
JJ = zeros(N*npde*npde*(d+1)*(d+1)+Nbf*npde*npde*d*d+len_b*len_b*npde*npde,1);
JAC_nz = zeros(N*npde*npde*(d+1)*(d+1)+Nbf*npde*npde*d*d+len_b*len_b*npde*npde,1);

ic = 0;
Indx = zeros(N,npde*(d+1));

UU = zeros(N,npde,d+1);
for j=1:(d+1)
   UU(:,:,j)= U((tri(:,j)),:);
end

for j1=1:(d+1)
for k1=1:npde
    
    h = (UU(:,k1,j1)+max(abs(UU(:,k1,j1)),1)*sqrt(eps))-UU(:,k1,j1);
    ss = (sign(UU(:,k1,j1))>= 0);
    h = (ss-(~ss)).*abs(h);
    
    UU(:,k1,j1) = UU(:,k1,j1) + h; 

    dU = zeros(N,d*npde);
    for i=1:npde
    for j=1:(d+1)
    for k=1:d
        dU(:,d*(i-1)+k) = dU(:,d*(i-1)+k) + Phi(:,(j-1)*d+k).*UU(:,i,j);
    end
    end
    end
    FF1 = zeros(N,npde*(d+1));
    for n = 1:nn
        xK = zeros(N,d);
        xdotK = zeros(N,d);
        uK = zeros(N,npde);
        udotK = zeros(N,npde);
        for j=1:(d+1)
            xK = xK + X(tri(:,j),:)*lambda(j,n);
            xdotK = xdotK + Xdot(tri(:,j),:)*lambda(j,n);
            uK = uK + UU(:,:,j)*lambda(j,n);
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
                FF1(:,npde*(j-1)+i) = FF1(:,npde*(j-1)+i) ...
                    + pdedef.volumeInt(dU,uK,udotK,dphiK,phiK,xK,t,i)*lambda(d+2,n);        
            end
        end
    end
    FF1 = bsxfun(@times,FF1,detE);  
    UU(:,k1,j1) = UU(:,k1,j1) - h; 
   
    h = 1.0./h;
    FF1 = FF1-FF;
    FF1 = bsxfun(@times,FF1,h);
    JAC_nz(ic+1:ic+N*npde*(d+1)) = reshape(FF1, N*npde*(d+1), 1);
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

% for boundary integrals (neumann bcs)

    % compute the areas/lengths of the boundary facets
    switch (d)
    case 1
        Ab = ones(Nbf,1);
    case 2
        Ab = sqrt(dot(X(tri_bf(:,2),:)-X(tri_bf(:,1),:), ...
                      X(tri_bf(:,2),:)-X(tri_bf(:,1),:),2));
    case 3
        Eb = cross(X(tri_bf(:,2),:)-X(tri_bf(:,1),:), ...
                   X(tri_bf(:,3),:)-X(tri_bf(:,1),:),2);
        Ab = sqrt(dot(Eb,Eb,2))*0.5;
    end
    % compute the integration weights
    [lambdab,nnb] = MovFEM_int_weights(d-1);
       
    % compute the residues
       
    FFb = zeros(Nbf,npde*d);
    dU_bf = dU(pdedef.bfElements,:);
    for n = 1:nnb
        xK = zeros(Nbf,d);
        uK = zeros(Nbf,npde);
        for j=1:d
            xK = xK + X(tri_bf(:,j),:)*lambdab(j,n);
            uK = uK + U(tri_bf(:,j),:)*lambdab(j,n);
        end
        for j=1:d
            phiK = repmat(lambdab(j,n),Nbf,1);
            for i=1:npde
                FFb(:,npde*(j-1)+i) = FFb(:,npde*(j-1)+i) ...
        + pdedef.boundaryInt(dU_bf,uK,phiK,xK,t,i,pdedef.bfMark)*lambdab(d+1,n);
            end
        end
    end
    FFb = bsxfun(@times,FFb,Ab);

Indx = zeros(Nbf,npde*d); 
UU = zeros(Nbf,npde,d);
for j=1:d
   UU(:,:,j)= U((tri_bf(:,j)),:);
end
UUK = zeros(Nbf,npde,d+1);
for j=1:d+1
   UUK(:,:,j)= U((tri(pdedef.bfElements(:),j)),:);
end
for j1=1:d
for k1=1:npde
    
    h = (UU(:,k1,j1)+max(abs(UU(:,k1,j1)),1)*sqrt(eps))-UU(:,k1,j1);
    ss = (sign(UU(:,k1,j1))>= 0);
    h = (ss-(~ss)).*abs(h);

    UU(:,k1,j1) = UU(:,k1,j1) + h;
    
    UUK(:,k1,j1) = UUK(:,k1,k1) + h; 
    dU_bf1 = zeros(Nbf,d*npde);
    for i=1:npde
    for j=1:(d+1)
    for k=1:d
        dU_bf1(:,d*(i-1)+k) = dU_bf1(:,d*(i-1)+k) ...
               + Phi(pdedef.bfElements,(j-1)*d+k).*UUK(:,i,j);
    end
    end
    end

    FFb1 = zeros(Nbf,npde*d);
    for n = 1:nnb
        xK = zeros(Nbf,d);
        uK = zeros(Nbf,npde);
        for j=1:d
            xK = xK + X(tri_bf(:,j),:)*lambdab(j,n);
            uK = uK + UU(:,:,j)*lambdab(j,n);
        end
        for j=1:d
            phiK = repmat(lambdab(j,n),Nbf,1);
            for i=1:npde
                FFb1(:,npde*(j-1)+i) = FFb1(:,npde*(j-1)+i) ...
        + pdedef.boundaryInt(dU_bf,uK,phiK,xK,t,i,pdedef.bfMark)*lambdab(d+1,n);
            end
        end
    end
    FFb1 = bsxfun(@times,FFb1,Ab);  
    UU(:,k1,j1) = UU(:,k1,j1) - h;
    UUK(:,k1,j1) = UUK(:,k1,j1) - h; 
   
    h = 1.0./h;
    FFb1 = FFb1-FFb;
    FFb1 = bsxfun(@times,FFb1,h);
    JAC_nz(ic+1:ic+Nbf*npde*d) = reshape(FFb1, Nbf*npde*d, 1);
    for j=1:(npde*d)
        Indx(:,j) = npde*(tri_bf(:,j1)-1)+k1;
    end
    JJ(ic+1:ic+Nbf*npde*d) = reshape(Indx, Nbf*npde*d, 1);
    for j=1:d
    for k=1:npde
        Indx(:,(j-1)*npde+k) = npde*(tri_bf(:,j)-1)+k;
    end
    end
    II(ic+1:ic+Nbf*npde*d) = reshape(Indx, Nbf*npde*d, 1);
    ic = ic + Nbf*npde*d;
end
end

% save to sparse matrix

    JAC = sparse(II(1:ic), JJ(1:ic), JAC_nz(1:ic), Nv*npde, Nv*npde);
    
% for dirichlet boundary conditions

    bfMark = repmat(pdedef.bfMark,1,d);
    JAC = JAC';
    for ipde=1:npde
        nodes_tmp = tri_bf(pdedef.bftype(:,ipde)==1,:);
        nodes_tmp = reshape(nodes_tmp,[],1);
        nodes = unique(nodes_tmp);
        JAC(:,npde*(nodes-1)+ipde) = 0;
    end
    JAC = JAC';
    
    % find nonzeros of JAC
    
    [II, JJ, JAC_nz] = find(JAC);
    ic = length(JAC_nz);
    
    % compute JAC for boundry points
    
    F = zeros(len_b,1);
    F1 = zeros(len_b,1);
    hh = zeros(len_b,1);
    sss= zeros(len_b,1);
    
    for ipde=1:npde
        bfMark_tmp = bfMark(pdedef.bftype(:,ipde)==1,:);
        nodes_tmp = tri_bf(pdedef.bftype(:,ipde)==1,:);
        bfMark_tmp = reshape(bfMark_tmp,[],1);
        nodes_tmp = reshape(nodes_tmp,[],1);
        [nodes, iia, ~] = unique(nodes_tmp);
        bfMark_tmp = bfMark_tmp(iia);
        len = length(nodes);
        if (len==0)
           continue;
        end
        F(1:len) = pdedef.dirichletRes(U(nodes,:),X(nodes,:),t,ipde,bfMark_tmp);
        for ipde1=1:npde
           
           hh(1:len) = (U(nodes,ipde1)+max(abs(U(nodes,ipde1)),1)*sqrt(eps))-U(nodes,ipde1);
           sss(1:len) = (sign(U(nodes,ipde1))>= 0);
           hh(1:len) = (sss(1:len)-(~sss(1:len))).*abs(hh(1:len));
          
           U(nodes,ipde1) = U(nodes,ipde1) + hh(1:len);
           F1(1:len) = pdedef.dirichletRes(U(nodes,:),X(nodes,:),t,ipde,bfMark_tmp);
           U(nodes,ipde1) = U(nodes,ipde1) - hh(1:len);
           F1(1:len) = (F1(1:len)-F(1:len))./hh(1:len);
           II(ic+1:ic+len) = npde*(nodes-1)+ipde;
           JJ(ic+1:ic+len) = npde*(nodes-1)+ipde1;
           JAC_nz(ic+1:ic+len) = F1(1:len);
           ic = ic + len;
        end
    end
    JAC = sparse(II(1:ic), JJ(1:ic), JAC_nz(1:ic), Nv*npde, Nv*npde);
    
% end of MovFEM_jac()
