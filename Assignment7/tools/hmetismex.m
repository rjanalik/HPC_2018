function metismex
% HMETISMEX Establish an interface between hMETIS and Matlab
% 
% hMetis
% [map,edgecut] = metismex('PartRecursive',A,nparts);
% [map,edgecut] = metismex('PartKway',A,nparts);
%
% Note that error checking is not done: make sure A is structurally
% symmetric or it will crash.
%
% See also METISDICE METISPART
