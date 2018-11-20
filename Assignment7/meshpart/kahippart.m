function [p1,p2] = kahippart(A,xy,randseed)
% KaHIPpart : Partition a graph using KaHIP default method.
%
% p = kahippart(A) returns a list of the vertices on one side of
%     a partition obtained by KaHIP 2.0 applied to the graph of A.
%     
% Optional arguments:
%   kahippart(A,xy) draws a picture of the partitioned graph,
%                   using the rows of xy as vertex coordinates.
%   [p1,p2] = kahippart(...)   also returns the list of vertices
%                              on the other side of the partition.
%
% See also kahipmex.m 
%
% Dimosthenis Pasadakis & Drosos Kourounis 29 Nov 17

if nargin < 2 || isempty(xy)
    xy = 0;
end;
picture = max(size(xy)) > 1;

if nargin<3, randseed=-1; end

map = kahip('KaHIP',A,2);
[p1,p2] = other(map);

if picture
   gplotpart(A,xy,p1);
   title('KaHIP default Partition')
end;
