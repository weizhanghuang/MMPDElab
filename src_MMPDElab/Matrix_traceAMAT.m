function tr = Matrix_traceAMAT(A,M,m,n)
%
% usage: tr = Matrix_traceAMAT(A,M,m,n) or
%        tr = Matrix_traceAMAT(A,M) for square matrices A and M.
%
% this function computes trace(A*M*A') for m-by-n matrix A and n-by-n matrix B.
%
% A:    a matrix A = [A_11, ..., A_m1, ..., A_1n, ..., A_mn], of size Nv-by-m*n.
%       each row of A is an m-by-n matrix.
% M:    a matrix M = [M_11, ..., M_n1, ..., M_1n, ..., M_nn], of size Nv-by-n*n.
%       each row of B is an n-by-n matrix.
% tr:   (output) tr = trace(A*M*A'), of size Nv-by-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

if (~exist('m','var') || isempty(m) || ~exist('n','var') || isempty(n))
      m = size(A,2);
      m = sqrt(m);
      n = m;
end

C = Matrix_mult(A,M,m,n,n);
C = Matrix_mult(C,Matrix_AT(A),m,n,m);
tr = Matrix_trace(C);

% end of Matrix_traceAMAT()
