function C = Matrix_AT(A,m,n)
%
% usage: C = Matrix_AT(A,m,n) or
%        C = Matrix_AT(A) for square matrices.
%
% this function computes C = A'.
% A:    a matrix A = [A_11, ..., A_m1, ..., A_1n, ..., A_mn], of size Nv-by-m*n.
%       each row of A is an m-by-n matrix.
% C:    (output) C = [C_11, ..., C_n1, ..., C_1m, ..., C_nm], of size Nv-by-m*n.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

C = zeros(size(A));

if (~exist('m','var') || isempty(m) || ~exist('n','var') || isempty(n))
      m = size(A,2);
      m = sqrt(m);
      n = m;
end

for i = 1:m
for j = 1:n
   C(:,n*(i-1)+j) = A(:,m*(j-1)+i);
end
end

% end of Matrix_AT()
