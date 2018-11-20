function kahipmex
% KaHIPmex : Establish an interface between KaHIP and Matlab
% KaHIP v2.0
%
% The graph partitioning framework KaHIP -- Karlsruhe High Quality Partitioning.
%
% The graph partitioning problem asks for a division of a graph's node set into 
% k equally sized blocks such that the number of edges that run between the 
% blocks is minimized. KaHIP is a family of graph partitioning programs.
% This mex executable incorporates KaFFPa (Karlsruhe Fast Flow Partitioner),
% which is a multilevel graph partitioning algorithm.
%
%
% USAGE:
%
% [map,edgecut] = kahip('KaHIP_NE',A,nparts);
% [map,edgecut] = kahip('KaHIP',A,nparts);
%
% The output is optional.
%
% Error checking is not done: make sure that A is structurally
% symmetric or it will crash.
%
% See also kahippart.m, kahipnepart.m
%
% Dimosthenis Pasadakis & Drosos Kourounis 29 Nov 17
