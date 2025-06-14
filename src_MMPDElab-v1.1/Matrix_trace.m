function tr = Matrix_trace(A)
%
% usage: tr = Matrix_trace(A)
%
% this function computes the trace of square matrix A.
%
% A:    a square matrix A = [A_11, ..., A_d1, ..., A_1d, ..., A_dd],
%       of size Nv-by-d*d. each row of A is a d-by-d matrix.
% tr:   (output) tr = trace(A), of size Nv-by-1.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

[Nv,d] = size(A);
d = sqrt(d);

tr = zeros(Nv,1);
for i = 1:d
   tr = tr + A(:,d*(i-1)+i);
end

% end of Matrix_trace()
