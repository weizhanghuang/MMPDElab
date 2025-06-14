function [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M,Linf_norm,Xi_ref)
%
% usage: [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M)
%        [Qgeo,Qeq,Qali] = MovMesh_MeshQualMeasure(X,tri,M,Linf_norm,Xi_ref)
%
% this function computes the geometric, equidistribution, and alignment measures
% (in maximum norm or L2 norm in xi) for a mesh defined via (X, tri)
% according to the metric tensor.
%
% if Xref is defined, the quality of mesh X is evaluated against Xref.
% otherwise, the quality of mesh X is evaluated against equilateral simplexes.
%
% X:    the coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
% M:    the metric tensor defined at the vertices of X, of size Nv-by-d*d.
% Linf_norm: (optional) {true or false} defines which norm, L_inf or L_2,
%       is used to compute the quality measures. default: true.
% Xi_ref: (optional) the coordinates of vertices of the reference mesh,
%       of size Nv-by-d. default: equilateral simplex.
% Qgeo: (output) the geometric quality measure.
% Qeq:  (output) the equidistribution quality measure.
% Qali: (output) the alignment quality measure.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

    % set basic parameters

    if (~exist('Linf_norm','var') || isempty(Linf_norm))
        Linf_norm = true; 
    end
    if (~islogical(Linf_norm))
       Linf_norm = true; 
    end

    N = size(tri, 1);
    d = size(tri, 2) - 1;
    
    % compute the quality measures

    [E, ~, ~, ~] = Matrix_edge(X, tri);
    mK = Matrix_average(M,tri);
    detM = Matrix_det(mK);
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
       EcinvK = reshape(inv(EcK),1,[]);
       Ecinv = repmat(EcinvK,N,1);
    else
       [~,Ecinv,~,~] = Matrix_edge(Xi_ref, tri);
    end
    detEc = 1./Matrix_det(Ecinv);
    detEc = abs(detEc)/factorial(d);
    FK = Matrix_mult(E,Ecinv);
    FK = Matrix_AT(FK);
    detFK = abs(Matrix_det(FK));
    qgeo = sum(FK.*FK,2)./(d*detFK.^(2/d));
    qali = Matrix_traceAMAT(FK,mK)./(d*detFK.^(2/d).*detM.^(1/d));
    detM = detFK.*sqrt(detM);
    sigmah = sum(detM.*detEc)/sum(detEc);
    qeq = detM/sigmah;
    if (Linf_norm)
       Qgeo = max(qgeo);
       Qali = max(qali);
       Qeq = max(qeq);
    else
       detEcinv = abs(Matrix_det(Ecinv));
       Area_C = sum(detEcinv);
       Qgeo = sqrt(sum(qgeo.^2.*detEcinv)/Area_C);
       Qali = sqrt(sum(qali.^2.*detEcinv)/Area_C);
       Qeq = sqrt(sum(qeq.^2.*detEcinv)/Area_C);
    end

% end of MovMesh_MeshQualMeasure()
