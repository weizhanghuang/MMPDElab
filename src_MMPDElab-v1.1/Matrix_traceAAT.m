function tr = Matrix_traceAAT(A)
%
% usage: tr = Matrix_traceAAT(A)
%
% this function computes the trace(A*A').
%
% A:    matrix A = [A_11, ..., A_m1, ..., A_1n, ..., A_mn], of size Nv-by-m*n.
%       each row of A is an m-by-n matrix.
% tr:   (output) tr = trace(A*A'), of size Nv-by-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

tr = dot(A,A,2);

% end of Matrix_traceAAT()
