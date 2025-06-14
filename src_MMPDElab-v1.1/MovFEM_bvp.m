function Unew = MovFEM_bvp(U,X,tri,tri_bf,pdedef,nonlinearsolver,MaxIter,Tol)
%
% usage: Unew = MovFEM_bvp(U,X,tri,tri_bf,pdedef)
%        Unew = MovFEM_bvp(U,X,tri,tri_bf,pdedef,nonlinearsolver,MaxIter,Tol)
%
% this function solves BVPs of the form:
%
%   sum_(i=1:npde) int_Omega [ F_i(nabla(u), u, u_t, nabla(v_i), v_i, x, t) ] dx
%
%         + sum_(i=1:npde) int_(Gamma_N_i) [ G_i(nabla(u), u, v_i, x, t) ] ds  = 0,
%
%                                                    for all v_i in H_i (i=1:npde)
% 
%   Dirichlet BC: Res_i(u, x, t) = 0 on Gamma_D_i (i=1:npde)
%
% where u = (u_1, ..., u_npde) are the unknown variables, Gamma_D_i is the part of
% boundary where a Dirichlet BC is defined for u_i, and H_i is a subspace of
% H^1(Omega) with functions vanishing on Gamma_D_i. P1 conforming fem is used.
% t is not used in the computation.
%
% U:    U = (u_1, u_2, ..., u_npde) initial values for the unknown variables,
%       of size Nv-by-npde.
% X:    the coordinates of vertices of the current mesh, of size Nv-by-d.
% tri:  the connectivity for all meshes, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets for all meshes, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% pdedef: a structure used to define the PDE system in the weak form.
%       it has 5 fields for the definition.
%
%       pdedef.bfMark: of size Nbf-by-1, can be used to mark the boundary segments
%                      (boundary facets), which is useful to define boundary 
%                      conditions and trace the boundary.
%       pdedef.bftype: of size Nbf-by-npde, specifies types of boundary condition
%                      on boundary facets whose numbering is based on tri_bf.
%                      0: Neumann BCs and 1: Dirichlet BCs.
%                      for For example,
%                       pdedef.bftype(3,2) = 1: variable u_2 has a Dirichlet BC
%                                               on the 3rd boundary facet.
%                       pdedef.bftype(3,2) = 0: variable u_2 has a Neumann BC
%                                               on the 3rd boundary facet. 
%
%       pdedef.volumeInt(du, u, ut, dv, v, x, t, i): defines F_i in the weak form.
%
%       pdedef.boundaryInt(du, u, v, x, t, i, bfMark): defines G_i in the weak form.
%
%       pdedef.dirichletRes(u, x, t, i, bfMark): defines function Res_i.
%       
%
%   note: all field functions should be defined in vector operation form.
%           the dimensions of the arguments are:
%
%           du(:,1:npde*d), u(:,1:npde), ut(:,1:npde), dv(:,d), v(:), x(:)
%
%           du(:,d*(i-1)+k) = d u_i/d_x_k: i = 1:npde, k=1:d. namely,
%           du = [u^(1)_x, u^(1)_y, u^(1)_z, u^(2)_x, ...]
%
% nonlinearsolver: (optional input) method used for solving nonlinear
%       algebraic systems. choices include 'newtons', 'fsolve'.
%       defacult is 'fsolve'.
% MaxIter: (optional input) maximum iteration allowed. default: MaxIter = 300.
% Tol:  (optional input) tolerance for iteration. default: Tol = 1e-6.
% Unew: (output) new approximation for the unknown variables, of size Nv-by-npde.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

% set basic parameters

fprintf('\n--- begin MovFEM_bvp ---\n');

    if (nargin==5 || ~exist('nonlinearsolver','var') || isempty(nonlinearsolver))
       nonlinearsolver = 'fsolve';
    end
    nonlinearsolver = lower(nonlinearsolver);
    if (~strcmp(nonlinearsolver,'fsolve') && ~strcmp(nonlinearsolver,'newtons'))
       nonlinearsolver = 'fsolve';
    end
    if (~exist('MaxIter','var') || isempty(MaxIter))
       MaxIter = 300;
    end
    if (~exist('Tol','var') || isempty(Tol))
       Tol = 1e-6;
    end
   
    [Nv,d] = size(X);
    npde = size(U,2);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Nbf = size(tri_bf,1);
    if (~isfield(pdedef,'bfMark'))
       pdedef.bfMark = ones(Nbf,1);
    end
    % find the elements containing the boundary facets
    switch(d)
        case 1
            pdedef.bfElements = zeros(2,1);
            [pdedef.bfElements(1),~] = find(tri==tri_bf(1));
            [pdedef.bfElements(2),~] = find(tri==tri_bf(2));
        case 2
            TR = triangulation(tri, X);
            pdedef.bfElements = cell2mat(edgeAttachments(TR,tri_bf));
        case 3
            TR = triangulation(tri, X);
            Elem1 = edgeAttachments(TR,tri_bf(:,1),tri_bf(:,2));
            Elem2 = edgeAttachments(TR,tri_bf(:,2),tri_bf(:,3));
            Elem = cellfun(@intersect,Elem1,Elem2,'UniformOutput',false);
            Elem1 = edgeAttachments(TR,tri_bf(:,1),tri_bf(:,3));
            Elem2 = cellfun(@intersect,Elem1,Elem,'UniformOutput',false);
            pdedef.bfElements = cell2mat(Elem2);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
% solve the nonlinear systems

    u = reshape(U',Nv*npde,1);
    ut = zeros(Nv*npde,1);
    Xdot = zeros(Nv,d);
    tn = 0;
    
    TolFun = Tol;
    TolX = Tol;
    
switch (nonlinearsolver)

    case 'fsolve'

        disp('    -- using fsolve');    
        opts = optimset('MaxIter',MaxIter,'TolFun',TolFun,'TolX',TolX, ...
                       'Jacobian','on','Display','iter','DerivativeCheck','off');
        unew = fsolve(@(u)MovFEM_bvp_fsolve(u,ut,npde,X,Xdot,tri,tri_bf,tn,pdedef), ...
                       u,opts);
                       
    case 'newtons' % modified Newton's method

        disp('    -- using Newtons iteration');
        for n=1:MaxIter
            [F, A] = MovFEM_bvp_fsolve(u,ut,npde,X,Xdot,tri,tri_bf,tn,pdedef);
            %[F,flag] = bicgstabl(A,F,TolX*0.01,Mit);
            F = A\F;
            unew = u - F;
            u = unew;
            fprintf('    newtons iteration n = %d  %e %e\n', n, ...
                            norm(F)/sqrt(Nv*npde), norm(F)/(norm(u)+eps));
            if (norm(F)/sqrt(Nv*npde)<=TolFun), break; end
            if (norm(F)<=TolFun*norm(u)), break; end
         end
         if (n == MaxIter) % Newton iteration fails
            error('Newtons iteration fails. MaxIter = %d\n', MaxIter);
         end
end

    Unew = reshape(unew,npde,[])';
    
fprintf('--- end MovFEM_bvp ---\n');
                                
% end of MovFEM_bvp()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [F, JAC] = MovFEM_bvp_fsolve(u,ut,npde,X,Xdot,tri,tri_bf,tn,pdedef)

F = MovFEM_rhs(0,u,ut,npde,X,Xdot,tri,tri_bf,tn,pdedef);
JAC = MovFEM_jac(0,u,ut,npde,X,Xdot,tri,tri_bf,tn,pdedef);

% end of MovFEM_bvp_fsolve()
