function [XF,TriF,TriF_parent] = MovMesh_MeshUniformRefine(X,Tri,Level)
%
% usage: [XF,TriF,TriF_parent] = MovMesh_MeshUniformRefine(X,Tri,Level)
%
% this function refines a simplicial mesh uniformly (Level) levels.
% on each level, 1 element is refined into 2^d elements:
% 1D: 1 interval into 2
% 2D: 1 triangle into 4
% 3D: 1 tetrahedron into 8
%
% X:    coordinates of vertices of the mesh, of size Nv-by-d.
% Tri:  the connectivity of the mesh, of size N-by-(d+1).
% Level: positive integer.
% XF:   (output) coordinates of vertices of the final mesh, of size NvF-by-d.
% TriF: (output) the connectivity of the final mesh, of size NF-by-(d+1).
% TriF_parent: (output) shows which original element a new element is in,
%       of size NF-by-1. 
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

% basic parameters

   if (~exist('Level','var') || isempty(Level) || Level< 0)
      Level = 1;
   end

   [Nv,d] = size(X);
   N = size(Tri,1);
   
   if (Level==0)   
      XF = X;
      TriF = Tri;
      TriF_parent = (1:N)';
      return;   
   end
   
   edges = zeros(2^(d*Level)*N*d*(d+1)/2,2);
   P = zeros(2^(d*Level)*N*d*(d+1)/2,2);
   
   XF = zeros(Nv+N*d*(d+1)*(2^((Level+1)*d)-2^d)/2/(2^d-1),d);
   TriF = zeros(2^(d*Level)*N,d+1);
   TriF_parent = repmat((1:N)',2^(d*Level),1);
   
   XC = zeros(Nv+N*d*(d+1)*(2^((Level+1)*d)-2^d)/2/(2^d-1),d);
   TriC = zeros(2^(d*Level)*N,d+1);
   
% start computation
   
   XC = X;
   TriC = Tri;
   
for l=1:Level
   
   % compute the new vertices (without repeated ones)
   
   % collect all edges
   icount = 0;
   for i=1:d+1
   for j=i+1:d+1
      edges(icount+1:icount+N,:) = [TriC(:,i),TriC(:,j)];
      icount = icount + N;
   end
   end
   % remove repeated edges. P(i,:) is the numbering of the edges contained
   % in the i-th element. the new vertices are ordered in the same way as
   % the edges.
   edges = sort(edges(1:N*d*(d+1)/2,:),2);
   [edges,~,P] = unique(edges(1:N*d*(d+1)/2,:),'rows');
   P = reshape(P,N,d*(d+1)/2) + Nv;
   P = [TriC, P];
   N_new = size(edges,1);
   XF = [XC;0.5*(XC(edges(:,1),:)+XC(edges(:,2),:))];
      
   % compute the connectivity
   
   switch (d)
   
   case 1
   
      TriF(1:N,:) = [P(:,1),P(:,3)];
      TriF(N+1:2*N,:) = [P(:,3),P(:,2)];
      
   case 2
               
      TriF(1:N,:) = [P(:,1),P(:,4),P(:,5)];
      TriF(N+1:2*N,:) = [P(:,2),P(:,6),P(:,4)];
      TriF(2*N+1:3*N,:) = [P(:,3),P(:,5),P(:,6)];
      TriF(3*N+1:4*N,:) = [P(:,4),P(:,6),P(:,5)];
      
   case 3
   
      TriF(1:N,:) = [P(:,1),P(:,5),P(:,6),P(:,7)];
      TriF(N+1:2*N,:) = [P(:,2),P(:,5),P(:,8),P(:,9)];
      TriF(2*N+1:3*N,:) = [P(:,3),P(:,6),P(:,8),P(:,10)];
      TriF(3*N+1:4*N,:) = [P(:,4),P(:,7),P(:,9),P(:,10)];
      TriF(4*N+1:5*N,:) = [P(:,5),P(:,6),P(:,7),P(:,8)];
      TriF(5*N+1:6*N,:) = [P(:,6),P(:,7),P(:,8),P(:,10)];
      TriF(6*N+1:7*N,:) = [P(:,7),P(:,8),P(:,9),P(:,10)];
      TriF(7*N+1:8*N,:) = [P(:,5),P(:,7),P(:,8),P(:,9)];
   end
   
   if (l==Level) break; end
   
   XC = XF;
   TriC = TriF(1:2^d*N,:);
   N = size(TriC,1);
   Nv = size(XC,1);
end

% end of MovMesh_MeshUniformRefine()
