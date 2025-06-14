function sqrtA = Matrix_sqrtA(A)
%
% usage: sqrtA = Matrix_sqrtA(A)
%
% this function computes the square root of symmetric and positive definite
% matrix A, i.e., sqrt(A). the analytical formulas based on Cayley-Hamilton
% theorem (e.g., see Franca, Comput. Math. Appl. 18 (1989), 459-466) are used
% for d = 1, 2, and 3 and the eigen-decomposition is used for d > 3.
% 
% A:    a symmetric and positive definite matrix A = [A_11, ..., A_d1, ..., A_1d,
%       ..., A_dd], of size Nv-by-d*d. each row of A is a d-by-d matrix. 
% sqrtA: (output) sqrtA = sqrt(A), symmetric and positive definite matrix,
%       of size Nv-by-d*d.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

    [Nv,d] = size(A);
    d = sqrt(d);
    
    switch (d)
       case 1
          sqrtA = sqrt(A);
       case 2   
          IIU = sqrt(Matrix_det(A));    
          sqrtA = A;
          sqrtA(:,1) = sqrtA(:,1) + IIU;
          sqrtA(:,4) = sqrtA(:,4) + IIU;
          IU = sqrt(Matrix_trace(sqrtA));
          sqrtA = bsxfun(@rdivide,sqrtA,IU);
          sqrtA(any(~isfinite(sqrtA),2),:) = abs(A(any(~isfinite(sqrtA),2),:));
       case 3
          IIIU = sqrt(Matrix_det(A));
          sqrtA = A;
          IU = Matrix_trace(sqrtA);
          IIU = 0.5*(IU.*IU-dot(sqrtA,sqrtA,2));
          k = sqrt(abs(IU.*IU-3*IIU));
          lambda = acos((IU.*(2*IU.*IU-9*IIU)+27*IIIU.*IIIU)./(2*k.^3));
          lambda = sqrt(abs((IU+2*k.*cos(lambda/3))/3));
          IU = lambda+sqrt(IU-lambda.*lambda+2*IIIU./lambda);
          IIU = sqrt(abs(IIU+2*IU.*IIIU));
          sqrtA = -Matrix_mult(sqrtA,sqrtA)+bsxfun(@times,sqrtA,IU.*IU-IIU);
          sqrtA(:,1) = sqrtA(:,1)+IU.*IIIU;
          sqrtA(:,5) = sqrtA(:,5)+IU.*IIIU;
          sqrtA(:,9) = sqrtA(:,9)+IU.*IIIU;
          sqrtA = bsxfun(@rdivide,sqrtA,IU.*IIU-IIIU);
          sqrtA(any(~isfinite(sqrtA),2),:) = abs(A(any(~isfinite(sqrtA),2),:));
       otherwise
          sqrtA = zeros(Nv,d*d);
          for i=1:Nv
             AA = reshape(A(i,:),d,d);
             [V,D] = eig(AA);
             D = sqrt(D);
             sqrtA(i,:) = reshape(V*D*V',1,[]);
          end
   end   

% end of MovMesh_sqrtA()

