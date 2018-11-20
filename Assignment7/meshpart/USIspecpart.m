function [part1,part2,sep1,sep2] = USIspecpart(A,xy,ignore);
% SPECPART : Spectral partition of a graph.
%
% [part1,part2] = specpart(A) returns a partition of the n vertices
%                 of A into two lists part1 and part2 according to the
%                 spectral bisection algorithm of Simon et al:  
%                 Label the vertices with the components of the Fiedler vector
%                 (the second eigenvector of the Laplacian matrix) and partition
%                 them about the median value.
%
% [part1,part2,sep1,sep2] = specpart(.) also returns the separating edges.
%
% If vertex coordinates are given as a second argument,
% specpart(A,xy) draws a picture of the result.
%
% See also LAPLACIAN, FIEDLER.
%

if nargin < 2
    xy = 0;
end;
picture = 1;

disp(' ');
disp(' HPC 2018 course:   ');
disp(' Implement your own spectral partitioning  ');
disp(' ');

% <<<<This code is just a dummy implementation to generate a partitioning
n = size(A,1);
map = zeros(n,1);
map(1:round((n/2))) = 0; map((round((n/2))+1):n) = 1;
[p1,p2] = other(map);
gplotpart(A,xy,p1);
title('Spectral Partition (dummy) using the Fiedler Eigenvector');
% <<<<This code is just a dummy implementation to generate a partitioning
