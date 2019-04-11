function detA = Matrix_det(A)
%
% usage: detA = Matrix_det(A)
%
% this function computes the determinant of sqaure matrix A.
%
% A:    A = [A_11, ..., A_d1, ..., A_1d, ..., A_dd], of size Nv-by-d*d.
%       each row of A is a d-by-d matrix.
% detA: (output) detA = det(A), a vector of size Nv-by-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

[Nv,d] = size(A);
d = sqrt(d);

switch (d)
    case 1
        detA = A(:,1);
    case 2 % [a11, a21, a12, a22] 
        detA = A(:,1).*A(:,4)-A(:,2).*A(:,3);
    case 3 % [a11, a21, a31, a12, a22, a32, a13, a23, a33]
        detA = A(:,1).*(A(:,5).*A(:,9) - A(:,6).*A(:,8)) ...
            - A(:,2).*(A(:,4).*A(:,9) - A(:,6).*A(:,7)) ...
            + A(:,3).*(A(:,4).*A(:,8) - A(:,5).*A(:,7));
    otherwise
        detA = zeros(Nv,1);
        for i=1:Nv
             detA(i) = det(reshape(A(i,:),d,d));
        end
end

% end of Matrix_det()
