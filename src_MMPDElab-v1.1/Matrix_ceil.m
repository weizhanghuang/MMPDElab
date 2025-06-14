function M1 = Matrix_ceil(M, beta)
%
% usage: M1 = Matrix_ceil(M, beta)
%
% this function puts a ceiling on the eigenvalues of symmetric and positive
% definite matrix M: lambda_max(M) <= beta
%
% M:    a symmetric and positive definite matrix M = [M_11, ..., M_d1,
%       ..., M_1d, ..., M_dd], of size Nv-by-d*d. each row of M is
%       a d-by-d matrix.
% beta: a positive number.
% M1:   (output) a symmetric and positive definite matrix, of size Nv-by-d*d.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

    trM = Matrix_trace(M);
    trM = 1./sqrt(1+trM.*trM/(beta*beta));
    M1 = bsxfun(@times,M,trM);

% end of Matrix_ceil()
