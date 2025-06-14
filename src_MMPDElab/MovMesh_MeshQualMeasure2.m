function [Qmax,Ql2] = MovMesh_MeshQualMeasure2(X,tri,M,Xi_ref)
%
% usage: [Qmax,Ql2] = MovMesh_MeshQualMeasure2(X,tri,M)
%        [Qmax,Ql2] = MovMesh_MeshQualMeasure2(X,tri,M,Xi_ref)
%
% this function computes
%     Ql2 = sum_K |K| \| (sigma/|Omega_c)^{2/d} JMJ^T - I\|_F^2
%     Qmax = max_K \| (sigma/|Omega_c)^{2/d} JMJ^T - I\|_F^2.
%
% if Xi_ref is defined, the quality of mesh X is evaluated against Xi_ref.
% otherwise, the quality of mesh X is evaluated against equilateral simplexes.
%
% X:    the coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
% M:    the metric tensor defined at the vertices of X, of size Nv-by-d*d.
% Xi_ref: (optinal) coordinates of vertices of the reference mesh,
%       of size Nv-by-d.
% Qmax: (output) the maximum norm of the quality measure.
% Ql2:  (output) the L2 norm of the quality measure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

    % set basic parameters

    N = size(tri, 1);
    d = size(tri, 2) - 1;
    
    % compute the measure

    [~, Einv, detE, ~] = Matrix_edge(X, tri);
    detE = abs(detE)/factorial(d);
    mK = Matrix_average(M,tri);
    detM = Matrix_det(mK);
    Minv = Matrix_inv(mK,detM);
    if (~exist('Xi_ref','var') || isempty(Xi_ref))
       switch (d)
          case 1
             EcK = 1;
          case 2
             EcK = [1, 0.5; 0, sqrt(3)*0.5];
          case 3
             EcK = [0, -2, -2; -2, 0, -2; -2, -2, 0];
       end
       EcK = (factorial(d)/abs(det(EcK)))^(1/d)*EcK;
       EcK = EcK/N^(1/d);
       EcK = reshape(EcK,1,[]);
       Ec = repmat(EcK,N,1);
    else
       [Ec,~,~,~] = Matrix_edge(Xi_ref, tri);
    end
    detEc = Matrix_det(Ec);
    detEc = abs(detEc)/factorial(d);
    sigma = dot(detE,sqrt(detM))/(eps+sum(detEc));
    J = Matrix_mult(Ec,Einv);
    JMJ = Matrix_mult(J,Minv);
    JMJ = Matrix_mult(JMJ, Matrix_AT(J));
    JMJ = sigma^(2/d)*JMJ - repmat(reshape(eye(d,d),1,[]),N,1);
    detEc = dot(JMJ, JMJ, 2);
    Ql2 = sqrt(dot(detE, detEc));
    Qmax = sqrt(max(detEc));

% end of MovMesh_MeshQualMeasure2()
