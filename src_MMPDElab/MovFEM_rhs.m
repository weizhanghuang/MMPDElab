function F = MovFEM_rhs(t,u,udot,npde,Xn,Xdot,tri,tri_bf,tn,pdedef)
%
% usage: F = MovFEM_rhs(t,u,udot,npde,Xn,Xdot,tri,tri_bf,tn,pdedef)
%
% this function computes F(t,u,udot) for the system F(t,u,udot) = 0
% resulting from the P1 conforming fem discretization of IBVPs described
% in MovFEM().
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
% tri_bf: the the boundary facets for all meshes, with each row representing
%       a facet on the boundary and containing d vertex IDs,
%       of size Nbf-by-d.
% tn:   beginning time at the current step.
% pdedef: a structure used to define the PDE system in the weak form;
%       see MovFEM().
% F:    (output) rhs, of size Nv*npde-by-1, F(npde*(j-1)+k), k=1:npde, j=1:Nv.
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

% define the nodes and weights for numerical integration

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
   
% compute dU at elements: dU(:,(i-1)*d+(1:d)) = \nabla u_i

    dU = zeros(N,d*npde);
    for i=1:npde
    for j=1:(d+1)
    for k=1:d
        dU(:,d*(i-1)+k)=dU(:,d*(i-1)+k)+Phi(:,(j-1)*d+k).*U(tri(:,j),i);
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
   
    % assemble rhs at nodes

    F = zeros(Nv,npde);
    for j=1:(d+1)
    for i=1:npde
        v = accumarray([tri(:,j);Nv], [FF(:,npde*(j-1)+i);0]);
        F(:,i) = F(:,i) + v; 
    end
    end

% for boundary integrals (neumann bcs)

    dU_bf = dU(pdedef.bfElements,:);
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
   
    % assemble rhs at nodes

    for j=1:d
    for i=1:npde
        v = accumarray([tri_bf(:,j);Nv], [FFb(:,npde*(j-1)+i);0]);
        F(:,i) = F(:,i) + v; 
    end
    end
    
% for dirichlet boundary conditions

    bfMark = repmat(pdedef.bfMark,1,d);
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
        F(nodes,ipde)=pdedef.dirichletRes(U(nodes,:),X(nodes,:),t,ipde,bfMark_tmp);
    end
   
    F = reshape(F',Nv*npde,1);
   
% end of MovFEM_rhs()
