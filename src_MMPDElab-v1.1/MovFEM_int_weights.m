function [lambda,nn] = MovFEM_int_weights(d)
%
% usuage: [lambda,nn] = MovFEM_int_weights(d)
%
% this function defines the weights for numerical integration used
% in fem discretization.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

% define the nodes and weights for numerical integration

    switch (d)
    
    case 0
        nn = 1;
        lambda = [1;1];
    
    case 1
        nn = 2;
        a = 1/sqrt(3);
        lambda = [0.5*(1-a), 0.5*(1+a);
                  0.5*(1+a), 0.5*(1-a);
                  0.5, 0.5];
    case 2
        % 3 interior point-based
        nn = 3;
        a = 2.0/3.0;
        b = 1.0/6.0; 
        lambda = [a, b, b;
                  b, a, b;
                  b, b, a;
                  1.0/3.0, 1.0/3.0, 1.0/3.0];
    case 3
        nn = 4;
        a = (5+3*sqrt(5))/20.0;
        b = (5-sqrt(5))/20.0;
        lambda = [a, b, b, b;
                  b, a, b, b;
                  b, b, a, b;
                  b, b, b, a;
                  0.25, 0.25, 0.25, 0.25];
    end

% end of MovFEM_int_weights()
