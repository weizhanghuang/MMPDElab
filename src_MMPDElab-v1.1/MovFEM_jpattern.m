function JP = MovFEM_jpattern(tri,npde)
%
% usage: JP = MovFEM_jpattern(tri,npde)
%
% this function computes the Jacobian pattern based on mesh connectivity.
% it is assumed that unknown valriables are ordered together at each mesh point.
%
% tri:  the connectivity of the mesh, of size N-by-(d+1).
% npde: number of the components of the physical solution.
% JP:   (output) space matrix containing the pattern of the Jacobian matrix
%       of the p1 fem discretization of the pdes.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

    [N,d] = size(tri);
    d = d - 1;
    Nv = max(max(tri));
   switch (d)
      case 1
      JP = sparse(reshape(tri(:,1).',[],1), ...
                         reshape(tri(:,2).',[],1),1,Nv,Nv);
      case 2
      JP = sparse(reshape(tri(:,[1 1 2]).',[],1), ...
                         reshape(tri(:,[2 3 3]).',[],1),1,Nv,Nv);
      case 3
      JP = sparse(reshape(tri(:,[1 1 1 2 2 3]).',[],1), ...
                         reshape(tri(:,[2 3 4 3 4 4]).',[],1),1,Nv,Nv);
   end
   JP = JP + JP.' + speye(Nv,Nv);
   JP(JP>1) = 1;
   JP = kron(JP,ones(npde,npde));
    
% end of MovFEM_jpattern()
