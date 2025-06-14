function [X,tri] = MovMesh_rect2tri(x,y,job)
%
% usage: [X,tri] = MovMesh_rect2tri(x,y)
%        [X,tri] = MovMesh_rect2tri(x,y,job)
%
% this function converts a rectangular mesh (defined as a tensor product
% of x and y) into a triangular mesh.
%
% x:    points in x-coordinate, of size nx-by-1.
% y:    points in y-coordinate, of size ny-by-1.
% job:  job = 1: a rectangular cell into 4 triangular cells
%           = 2: a rectangular cell into 2 triangular cells (right)
%           = 3: a rectangular cell into 2 triangular cells (left)
% X:    (output) coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  (output) the connectivity of the mesh, of size N-by-(d+1).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

    if (~exist('job','var') || isempty(job))
        job = 1; 
    end
    if (job ~= 1 && job ~= 2 && job ~= 3)
       job = 1;
    end

   jmax = max(size(x));
   kmax = max(size(y));
   [xx, yy] = ndgrid(x,y); 
   X = [xx(:), yy(:)];
   N = (jmax-1)*(kmax-1);
   
   % nodes: 1(0,0), 2(1,0), 3(1,1), 4(0,1)
   % range of vertices: jmax*(k-1)+j
   
   I1 = (1:(jmax*(kmax-1)))';
   I1(jmax:jmax:(jmax*(kmax-1))) = [];

   I2 = I1 + 1;
   I3 = I1 + jmax + 1;
   I4 = I1 + jmax;

   switch (job)

   case 1 % into 4 triangles

      tri = zeros(4*N,3);
      % for middle points
      xm = 0.25*(X(I1,1)+X(I2,1)+X(I3,1)+X(I4,1));
      ym = 0.25*(X(I1,2)+X(I2,2)+X(I3,2)+X(I4,2));
      X = [X; xm(:), ym(:)];
      I5 = jmax*kmax + (1:N)';
      % for element 1-2-5
      tri(1:N,:) = [I1, I2, I5];
      % for element 2-3-5
      tri(N+1:2*N,:) = [I2, I3, I5];
      % for element 3-4-5
      tri(2*N+1:3*N,:) = [I3, I4, I5];
      % for element 1-4-5
      tri(3*N+1:4*N,:) = [I1, I5, I4];
      
   case 2 % into 2 triangles (right)

      tri = zeros(2*N,3);
      % for element 1-2-3
      tri(1:N,:) = [I1, I2, I3];
      % for element 1-3-4
      tri(N+1:2*N,:) = [I1, I3, I4];

   case 3 % into 2 triangles (left)

      tri = zeros(2*N,3);
      % for element 1-2-4
      tri(1:N,:) = [I1, I2, I4];
      % for element 2-3-4
      tri(N+1:2*N,:) = [I2, I3, I4];
   end

% end of MovMesh_rect2tri()
