function fnew = MovMesh_LinInterp(f,X,QP,tri,tri_bf,useDelaunayTri)
%
% usage: fnew = MovMesh_LinInterp(f,X,QP,tri,tri_bf)
%        fnew = MovMesh_LinInterp(f,X,QP,tri,tri_bf,useDelaunayTri)
%
% this function performs linear interpolation of f (defined on X)
% at query points QP using triangulation or Delaunay triangulation.
%
% Note: this funtion may not may properly with non-convex domains, especially
%       when useDelaunayTri = true (with delaunayTriangulation clas).
%       The pointLocation function of triangulation class is much SLOWER than
%       the pointLocation function of delaunayTriangulation class in general.
%       However, the latter may give undesired results when the domain is
%       not convex since delaunayTriangulation is typically for convex domains
%       only.
%       
% f:    function values defined at vertices of the mesh (X,tri),
%       of size Nv-by-Ncomp, with Ncomp being the number of componoents.
%       f(i,:) = f(x_i).
% X:    the coordinates of vertices of the mesh, of size Nv-by-d.
% QP:   query points, of size Nq-by-d, with Nq being the number of query points.
% tri:  the connectivity of the mesh, of size N-by-(d+1).
% tri_bf: the boundary facets of the mesh, with each row representing
%       a facet on the boundary and containing d vertex IDs, of size Nbf-by-d.
% useDelaunayTri: {true, false}.
%       default is true, i.e. use delaunayTriangulation class.
%       otherwise, use triangulation class.
% fnew: (output) function values at QP, of size Nq-by-Ncomp.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2019  Weizhang Huang (whuang@ku.edu)

   % set basic parameters

   if (~exist('useDelaunayTri','var') || isempty(useDelaunayTri))
      useDelaunayTri = true;
   end
   if (~islogical(useDelaunayTri))
      useDelaunayTri = true;
   end
   Ncomp = size(f,2);
   [Nq,d] = size(QP);
   fnew = zeros(Nq,Ncomp);

   % perform linear interpolation when d = 1:

   if (d==1)
      fnew = interp1(X,f,QP,'linear');
      return;
   end

   % form a delaunay triangulation/or a general triangulation
   % Delaunay triangulation is fast, good for convex domains.
   % For non-convex domain, this option may lead to boundary points moving
   % out of domain when used as linear interpolation for obtaining
   % the new physical mesh from the new computational mesh.

   if (useDelaunayTri) 
      TR = delaunayTriangulation(X);
   else 
      TR = triangulation(tri,X);
   end

   % find the enclosing triangles or tetrahedra and barycentric coordinates
   % for all query points in QP

   [ti,B] = pointLocation(TR,QP);

   % check and treat the situation when some points lie outside of convex hull

   nodes_o = find(isnan(ti));
   N_o = length(nodes_o);

   if (N_o>0)
      vi = nearestNeighbor(TR,QP(nodes_o,:)); % find closest vertices
      TI = vertexAttachments(TR,vi);
      
      QQP = mat2cell(QP(nodes_o,:),ones(1,N_o),d);
      [groupMat,groupj] = cellfun(@(A1,A2)find_bc(A1,A2,TR),TI,QQP, ...
                                  'UniformOutput',false);
      B(nodes_o,:) = cell2mat(groupMat);
      ti(nodes_o,:) = cell2mat(groupj);
      %{ 
      for i=1:N_o
         len = length(TI{i});
         BB = cartesianToBarycentric(TR,(TI{i})',repmat(QP(nodes_o(i),:),len,1));
         % [~,j] = max(min(BB,[],2));
         [~,j] = min(abs(sum(BB,2)-sum(abs(BB),2)));
         B(nodes_o(i),:) = BB(j,:);
         ti(nodes_o(i),:) = TI{i}(j);
      end
      %}
   end

   % compute linear interpolation

   nodes = TR.ConnectivityList(ti,:);
   for j=1:Ncomp
      ff = f(:,j);
      triF = ff(nodes);
      fnew(:,j) = dot(B,triF,2);
   end

% end of MovMesh_LinInterp()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,j] = find_bc(TI,QP,TR)

    len = length(TI(:));
    BB = cartesianToBarycentric(TR,TI(:),repmat(QP(:)',len,1));
    % [~,j] = max(min(BB,[],2));
    [~,j] = min(abs(sum(BB,2)-sum(abs(BB),2)));
    B = BB(j,:);
    j = TI(j);

% end of find_bc()
