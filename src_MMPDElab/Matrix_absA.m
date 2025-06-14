function absA = Matrix_absA(A)
%
% usage: absA = Matrix_absA(A)
%
% this function computes the absolute of square matrix A, i.e., sqrt(A*A).
% the analytical formulas based on Cayley-Hamilton theorem
% (e.g., see Franca, Comput. Math. Appl. 18 (1989), 459-466) are used
% for d = 1, 2, and 3 and the eigen-decomposition is used for d > 3.
% 
% A:    a symmetric square matrix A = [A_11, ..., A_d1, ..., A_1d, ..., A_dd],
%       of size Nv-by-d*d. each row of A is a d-by-d matrix.
% absA: (output) absA = sqrt(A*A), symmetric matrix, of size Nv-by-d*d.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)
%
% June 14, 2021: the analytical formula given in Franca (1989) does not seem
% robust for certain matrices for d=3. Moreover, it may give a square root
% that may not be wanted. for this reason, we switch back to
% eigen-decomposition for d=3.
%

    [Nv,d] = size(A);
    d = sqrt(d);
    
    switch (d)
    
       case 1
          absA = abs(A);
       case 2
          IIU = abs(Matrix_det(A));
          absA = Matrix_mult(A,A);
          absA(:,1) = absA(:,1) + IIU;
          absA(:,4) = absA(:,4) + IIU;
          IU = sqrt(Matrix_trace(absA));
          absA = bsxfun(@rdivide,absA,IU);
          absA(any(~isfinite(absA),2),:) = abs(A(any(~isfinite(absA),2),:));
%{
       case 3
          IIIU = abs(Matrix_det(A));
          absA = Matrix_mult(A,A);
          IU = Matrix_trace(absA);
          IIU = 0.5*(IU.*IU-dot(absA,absA,2));
          k = sqrt(abs(IU.*IU-3*IIU));
          lambda = acos((IU.*(2*IU.*IU-9*IIU)+27*IIIU.*IIIU)./(2*k.^3));
          lambda = sqrt(abs((IU+2*k.*cos(lambda/3))/3));
          IU = lambda+sqrt(IU-lambda.*lambda+2*IIIU./lambda);
          IIU = sqrt(abs(IIU+2*IU.*IIIU));
          absA = -Matrix_mult(absA,absA)+bsxfun(@times,absA,IU.*IU-IIU);
          absA(:,1) = absA(:,1)+IU.*IIIU;
          absA(:,5) = absA(:,5)+IU.*IIIU;
          absA(:,9) = absA(:,9)+IU.*IIIU;
          absA = bsxfun(@rdivide,absA,IU.*IIU-IIIU);
          absA(any(~isfinite(absA),2),:) = abs(A(any(~isfinite(absA),2),:));
%}
       otherwise
          absA = zeros(Nv,d*d);
          for i=1:Nv
             [V,D] = eig(reshape(A(i,:),d,d));
             absA(i,:) = reshape(V*abs(D)*V',1,[]);
          end
   end    

% end of Matrix_absA()
