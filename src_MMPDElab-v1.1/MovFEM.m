function [Unew,dt0,dt1] = MovFEM(t,dt,U,X,Xdot,tri,tri_bf,pdedef, ...
                      fixed_step,relTol,absTol,direct_ls,ControlWeights)
%
% usage: [Unew,dt0,dt1] = MovFEM(t,dt,U,X,Xdot,tri,tri_bf,pdedef)
%        [Unew,dt0,dt1] = MovFEM(t,dt,U,X,Xdot,tri,tri_bf,pdedef, ...
%                            fixed_step,relTol,absTol,direct_ls,ControlWeights)
%
% this function solves IBVPs of the form:
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
% the resulting ode system is integrated using the fifth-order Radau IIA method.
%
% t:    current time.
% dt:   intended time step for integrating physical pdes.
% U:    U = (u_1, u_2, ..., u_npde) current solution, of size Nv-by-npde.
% X:    the coordinates of vertices of the current mesh, of size Nv-by-d.
% Xdot: mesh velocity, of size Nv-by-d.
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
% fixed_step: (optional input) true or false. default is false.
% relTol, absTol: (optional input) relative and absolute tolerances for
%       time step selection. default are relTol = 1e-4, absTol = 1e-6.
% direct_ls: (optional input) true or false. default is true (direct sparse
%       solver for linear systems is used).
% ControlWeights: (optional input) nonnegative vector of size (Nv*npde)-by-1.
%       it defines the weights of the components of the solution used in the error
%       estimation for time step selection.
% Unew: (output) new approximation for the unknown variables, of size Nv-by-npde.
% dt0:  (output) time step size used in the last step in integrating physical pdes.
% dt1:  (output) time step size predicted for the next step.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

% set basic parameters

fprintf('\n--- begin MovFEM ---\n');

    [Nv,d] = size(X);
    npde = size(U,2);
    u = reshape(U',Nv*npde,1);
    tn = t;

    % check variables
    
    if (~exist('fixed_step','var') || isempty(fixed_step))
        fixed_step = false;
        relTol = 1e-4;
        absTol = 1e-6;
    end
    if (~exist('relTol','var') || isempty(relTol))
        relTol = 1e-4;
    end
    if (~exist('absTol','var') || isempty(absTol))
        absTol = 1e-6;
    end
    if (~exist('direct_ls','var') || isempty(direct_ls) || direct_ls)
        direct_ls = true;
    end
    if (~islogical(direct_ls))
       direct_ls = true; 
    end
    if (~exist('ControlWeights','var') || isempty(ControlWeights) ...
                                  || length(ControlWeights)<(npde*Nv))
        ControlWeights = ones(npde*Nv,1);
    end
    
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

    
    % set initial approxination for ut
    ut = zeros(Nv*npde,1);   
    
    disp('    -- using RADAU45TS --');
    [unew,dt0,dt1] = MovFEM_radau45TS(@MovFEM_rhs,@MovFEM_implicit_jac,t,dt,u,ut, ...
                           relTol,absTol,fixed_step,direct_ls,ControlWeights, ...
                           npde,X,Xdot,tri,tri_bf,tn,pdedef);
    Unew = reshape(unew,npde,Nv)';
    
fprintf('--- end of MovFEM ---\n');
                                
% end of MovFEM()
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dfdy,dfdyt] = MovFEM_implicit_jac(t,u,ut,npde,X,Xdot,tri,tri_bf,tn,pdedef)
  
   dfdy = MovFEM_jac(t, u, ut, npde, X, Xdot, tri, tri_bf, tn, pdedef);
   dfdyt = MovFEM_mass(t, u, ut, npde, X, Xdot, tri, tri_bf, tn, pdedef);
   
% end of MovFEM_implicit_jac()
