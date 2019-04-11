function C = Matrix_mult(A,B,m,n,l)
%
% usage: C = Matrix_mult(A,B,m,n,l) or
%        C = Matrix_mult(A,B) for square matrices A and B.
%
% this function computes A*B, where m-by-n matrix A and n-by-l matrix B.
%
% A:    a matrix A = [A_11, ..., A_m1, ..., A_1n, ..., A_mn], of size Nv-by-m*n.
%       each row of A is an m-by-n matrix.
% B:    a matrix B = [B_11, ..., B_n1, ..., B_1l, ..., B_nl], of size Nv-by-n*l.
%       each row of B is an n-by-l matrix.
% C:    (output) C = [C_11, ..., C_m1, ..., C_1l, ..., C_ml], of size Nv-by-m*l.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

if (~exist('m','var') || isempty(m) || ~exist('n','var') || isempty(n) ...
    || ~exist('l','var') || isempty(l))
      m = size(A,2);
      m = sqrt(m);
      n = m;
      l = m;
end

C = zeros(size(A,1),m*l);

for i=1:m
for k=1:l
for j=1:n
   C(:,m*(k-1)+i) = C(:,m*(k-1)+i)+A(:,m*(j-1)+i).*B(:,n*(k-1)+j);
end
end
end

% end of Matrix_mult()
