function Ainv = Matrix_inv(A,detA)
%
% usage: Ainv = Matrix_inv(A,detA)
%
% this function computes the inverse of sqaure matrix A.
%
% A:    A = [A_11, ..., A_d1, ..., A_1d, ..., A_dd], of size Nv-by-d*d.
%       each row of A is a d-by-d matrix.
% Ainv: (output) Ainv = inv(A), square matrix of size Nv-by-d*d.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

if (nargin==1)
   detA = Matrix_det(A);
end

[Nv,d] = size(A);
d = sqrt(d);

switch (d)
    case 1
        Ainv(:,1) = 1./A(:,1);
    case 2 % [a11, a21, a12, a22]
        Ainv(:,1) =  A(:,4)./detA(:);
        Ainv(:,2) = -A(:,2)./detA(:);
        Ainv(:,3) = -A(:,3)./detA(:);
        Ainv(:,4) =  A(:,1)./detA(:);
    case 3 % [a11, a21, a31, a12, a22, a32, a13, a23, a33]
        Ainv(:,1) = (A(:,5).*A(:,9) - A(:,6).*A(:,8))./detA(:);
        Ainv(:,2) = (A(:,3).*A(:,8) - A(:,2).*A(:,9))./detA(:);
        Ainv(:,3) = (A(:,2).*A(:,6) - A(:,3).*A(:,5))./detA(:);
        Ainv(:,4) = (A(:,6).*A(:,7) - A(:,4).*A(:,9))./detA(:);
        Ainv(:,5) = (A(:,1).*A(:,9) - A(:,3).*A(:,7))./detA(:);
        Ainv(:,6) = (A(:,3).*A(:,4) - A(:,1).*A(:,6))./detA(:);
        Ainv(:,7) = (A(:,4).*A(:,8) - A(:,5).*A(:,7))./detA(:);
        Ainv(:,8) = (A(:,2).*A(:,7) - A(:,1).*A(:,8))./detA(:);
        Ainv(:,9) = (A(:,1).*A(:,5) - A(:,2).*A(:,4))./detA(:);
    otherwise
        Ainv = zeros(Nv,d*d);
        for i=1:Nv
             AA = reshape(A(i,:),d,d);
             Ainv(i,:) = reshape(inv(AA),d,d);
        end
end

% end of Matrix_inv()
