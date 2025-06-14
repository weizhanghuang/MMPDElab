function err = MovFEM_Error_P1Linf(uexact,t,X,U,tri,tri_bf,varargin)
%
% usage: err = MovFEM_Error_P1Linf(uexact,t,X,U,tri,tri_bf,varargin)
%
% this function computes the L-inf error in P1 finite element approximation.
%
% uexact: defines the exact solution of form U = uexact(t,x,varargin),
%       output U should be of size Nv-by-npde.
% t:    current time.
% X:    coordinates of mesh vertices, of size Nv-by-d.
% U:    U = (u_1, u_2, ..., u_npde) the unknown variables, of size Nv-by-npde.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
%       tri(i,:) contains IDs of all vertices in element i.
% tri_bf: the boundary facets of the mesh, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

% compute basic parameters

    d = size(X,2);

% define integration points

    nn = 6;
    lambda = zeros(nn,d);
    for i=1:d
       lambda(:,i) = linspace(0,1,nn)';
    end
   
% compute the error

   err = 0.0;

   switch (d)
   
   case 1   
      ic = 0;
      for i=1:nn
         lambda1 = 1-lambda(i,1);
         if lambda1 < 0, continue; end
         ic = ic+1;
         xx = lambda(i,1)*X(tri(:,1),:) ...
            + lambda1    *X(tri(:,2),:);
         UU = lambda(i,1)*U(tri(:,1),:) ...
            + lambda1    *U(tri(:,2),:);
         Ue = uexact(t, xx, varargin{:});
         Ue = Ue-UU;
         err = max(norm(Ue,Inf),err);
      end
   case 2
      ic = 0;
      for i=1:nn
      for j=1:nn
         lambda1 = 1-(lambda(i,1)+lambda(j,2));
         if lambda1 < 0, continue; end
         ic = ic+1;
         xx = lambda(i,1)*X(tri(:,1),:) ...
            + lambda(j,2)*X(tri(:,2),:) ...
            + lambda1    *X(tri(:,3),:);
         UU = lambda(i,1)*U(tri(:,1),:) ...
            + lambda(j,2)*U(tri(:,2),:) ...
            + lambda1    *U(tri(:,3),:);
         Ue = uexact(t, xx, varargin{:});
         Ue = Ue-UU;
         err = max(norm(Ue,Inf),err);
      end
      end   
   case 3
      ic = 0;
      for i=1:nn
      for j=1:nn
      for k=1:nn
         lambda1 = 1-(lambda(i,1)+lambda(j,2)+lambda(k,3));
         if lambda1 < 0, continue; end
         ic = ic+1;
         xx = lambda(i,1)*X(tri(:,1),:) ...
            + lambda(j,2)*X(tri(:,2),:) ...
            + lambda(k,3)*X(tri(:,3),:) ...
            + lambda1    *X(tri(:,4),:);
         UU = lambda(i,1)*U(tri(:,1),:) ...
            + lambda(j,2)*U(tri(:,2),:) ...
            + lambda(k,3)*U(tri(:,3),:) ...
            + lambda1    *U(tri(:,4),:);
         Ue = uexact(t, xx, varargin{:});
         Ue = Ue-UU;
         err = max(norm(Ue,Inf),err);
      end
      end
      end     
   end
      
% end of MovFEM_Error_P1Linf()
