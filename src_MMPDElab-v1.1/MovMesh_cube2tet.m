function [X,tri] = MovMesh_cube2tet(x,y,z)
%
% usage: [X,tri] = MovMesh_cube2tet(x,y,z)
%
% this function converts a cuboid mesh (defined as a tensor product
% by x and y and z) into a tetrahedral mesh (each subcuboid is divided
% into 6 tetrahedra).
%
% x:    points in x-coordinate, of size nx-by-1.
% y:    points in y-coordinate, of size ny-by-1.
% z:    points in z-coordinate, of size nz-by-1.
% X:    (output) coordinates of vertices of the mesh, of size Nv-by-d.
% tri:  (output) the connectivity of the mesh, of size N-by-(d+1).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   imax = max(size(x));
   jmax = max(size(y));
   kmax = max(size(z));

   N = (imax-1)*(jmax-1)*(kmax-1);  
   [xx, yy, zz] = ndgrid(x,y,z); 
   X = [xx(:), yy(:), zz(:)];
   tri = zeros(6*N,4);
   
   % nodes: 1(0,0,0), 2(0,0,1), 3(1,0,0), 4(1,0,1),
   %        5(0,1,0), 6(0,1,1), 7(1,1,0), 8(1,1,1)
   % range of vertices: imax*( jmax*(k-1)+(j-1) )+i

   I1 = find(X(:,1)<x(end)-10*eps&X(:,2)<y(end)-10*eps&X(:,3)<z(end)-10*eps);
   I2 = I1 + imax*jmax;
   I3 = I1 + 1;
   I4 = I1 + imax*jmax + 1;
   I5 = I1 + imax;
   I6 = I1 + imax*jmax + imax;
   I7 = I1 + imax + 1;
   I8 = I1 + imax*jmax + imax + 1;

   % Element 1: 1,8,3,7
   tri(1:N,:) = [I1, I8, I3, I7];
   % Element 2: 1,4,3,8
   tri(N+1:2*N,:) = [I1, I4, I3, I8];
   % Element 3: 1,8,2,4
   tri(2*N+1:3*N,:) = [I1, I8, I2, I4];   
   % Element 4: 1,8,7,5
   tri(3*N+1:4*N,:) = [I1, I8, I7, I5];   
   % Element 5: 1,8,5,6
   tri(4*N+1:5*N,:) = [I1, I8, I5, I6];   
   % Element 6: 1,6,2,8
   tri(5*N+1:6*N,:) = [I1, I6, I2, I8];

% end of MovMesh_cube2tet()
