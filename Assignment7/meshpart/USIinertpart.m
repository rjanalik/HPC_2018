function p = USIinertpart(A,xy,picture);
% INERTPART : Inertial partition of a graph.
%
% p = USIinertpart(A,xy) returns a list of the vertices on one side of a partition
%     obtained by bisection with a line to a moment of inertia
%     of the vertices, considered as points in Euclidean space.
%     Input A is the adjacency matrix of the mesh (used only for the picture!);
%     each row of xy is the coordinates of a point in d-space.
%
% USIinertpart(A,xy,1) also draws a picture.
%

if nargin < 3
    picture = (nargout == 0);
end;
[n,d] = size(xy);

disp(' ');
disp(' HPC 2018 course:   ');
disp(' Implement your own inertial bisection partitioning  ');
disp(' ');

% <<<<This code is just a dummy implementation to generate a partitioning
n = size(A,1);
map = zeros(n,1);
map(1:round((n/2))) = 0; map((round((n/2))+1):n) = 1;
[p1,p2] = other(map);
gplotpart(A,xy,p1);
title('Inertial Partition (dummy)');
disp(' Here we will generate a dummy partitioning ...');
% <<<<This code is just a dummy implementation to generate a partitioning
