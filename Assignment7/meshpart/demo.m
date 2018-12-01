function meshdemo(whichdemos)
% MESHDEMO : Demo of mesh partitioning toolkit.
%
% meshdemo(whichdemos);
% The argument is a demo or a list of demos to run.  
% The default is to run them all.

% John Gilbert, Xerox PARC, August 1994.
% Copyright (c) 1990-1996 by Xerox Corporation.  All rights reserved.
% HELP COPYRIGHT for complete copyright and licensing notice.
%
% Modified by Olaf Schenk, for HPC class at USI Lugano

addpath('../tools');

clear all
close all

warning('off','all');
warning

if nargin < 1
    whichdemos = [1 2 3];
end;

format compact;

disp('          *********************************************')
disp('          *** Mesh Partitioning and Separator Demos ***');
disp('          *********************************************')
disp(' ');
disp(' The file "meshes.mat" contains some sample meshes with coordinates.');
disp(' ');
load meshes;
whos;


for nmesh = 1:9
close all; clf reset;

if (nmesh==1)  
   disp(' ');
   disp(' Function "grid5rec" produces a rectangular grid:');
   disp(' ');
   disp('[A,xy] = grid5rec(8,80);');
   disp(' ');
   [A,xy] = grid5rec(8, 80);
   Tmesh = A;
   Tmeshxy= xy;
end
if (nmesh==2)  
   disp(' ');
   disp(' Function "grid5rec" produces a rectangular grid:');
   disp(' ');
   disp('[A,xy] = grid5rec(80,8);');
   disp(' ');
   [A,xy] = grid5rec(80, 8);
   Tmesh = A;
   Tmeshxy= xy;
end

if (nmesh==3)  
   disp(' ');
   disp(' Function "gridt" produces a triangular grid:');
   disp(' ');
   disp(' (See also grid5, grid7, grid9, grid3d, grid3dt.)');
   disp(' ');
   disp('[A,xy] = gridt(20);');
   disp(' ');
   [A,xy] = gridt(20);
   Tmesh = A;
   Tmeshxy= xy;
end
if (nmesh==4)  
   disp(' ');
   disp(' Function "gridt" produces a triangular grid:');
   disp(' ');
   disp(' (See also grid5, grid7, grid9, grid3d, grid3dt.)');
   disp(' ');
   disp('[A,xy] = gridt(20);');
   disp(' ');
   [A,xy] = grid9(30);
   Tmesh = A;
   Tmeshxy= xy;
end
if (nmesh==5)  
   Tmesh = Smallmesh;
   Tmeshxy= Sxy;
end
if (nmesh==6)  
   Tmesh = Tapir;
   Tmeshxy= Txy;
end
if (nmesh==7)  
   Tmesh = Eppstein; 
   Tmeshxy= Exy;
end
if (nmesh==8)  
   disp(' ');
   disp(' This demo shows the finite element mesh for a NASA airfoil   ');
   disp(' ');
   Tmesh = Airfoil; 
   Tmeshxy= Airfoilxy;
end
if (nmesh==9)  
   disp(' ');
   disp(' Generate a mesh (with 6*k points) whose best edge separator ');
   disp(' has size 2, but for which the spectral algorithm gives a    ');
   disp(' separator of size O(k). From Guattery and Miller, "On the   ');
   disp(' performance of spectral graph partitioning methods,"        ');
   disp(' SODA 1995.)                                                 ');
   disp(' ');
   disp('[A,xy = cockroach(60);');
   disp(' ');
   [A,xy] = cockroach(60);
   Tmesh = A;
   Tmeshxy= xy;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(whichdemos==1) % Various partitioning methods on all meshes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('          *********************************************')
disp('          ***      Various Partitioning Methods    *** ');
disp('          *********************************************')
disp(' ');
disp(' ');
    

if (nmesh==1)  
   disp('An initial rectangular  grid5rec(8,80) mesh');
end
if (nmesh==2)  
   disp('An initial rectangular  grid5rec(80,8) mesh');
end
if (nmesh==3)  
   disp('  gridt(20) mesh');
end
if (nmesh==4)  
   disp('  gridt9(20) mesh');
end
if (nmesh==5)  
   disp(' Small mesh ');
end
if (nmesh==6)  
   disp(' "Tapir" is a test of a no-obtuse-angles mesh generation algorithm');
   disp(' due to Bern, Mitchell, and Ruppert.  ');
end
if (nmesh==7)  
   disp('  Eppstein mesh');
end
if (nmesh==8)  
   disp(' ');
   disp('A mesh around the outside of an airfoil                           ');
   disp(' ');
end
if (nmesh==9)  
   disp(' ');
   disp(' A bad mesk for the spectral algorithm                       ');
   disp(' ');
end


figure(1)
disp('gplotg(Tmesh,Tmeshxy);');
disp('  ');
gplotg(Tmesh,Tmeshxy);
nvtx = size(Tmesh,1);
nedge = (nnz(Tmesh)-nvtx)/2;
xlabel([int2str(nvtx) ' vertices, ' int2str(nedge) ' edges'],'visible','on');

disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' We''ll try a few different partitioners on the mesh:');
disp(' ');
disp(' ');
disp(' First is Coordinate bisection partition of a mesh.  ');
disp(' p = coordpart(A,xy) returns a list of the vertices  ');
disp(' on one side of a partition obtained by bisection    ');
disp(' perpendicular to a coordinate axis.  We try every   ');
disp(' coordinate axis and return the best cut. Input A is '); 
disp(' the adjacency matrix of the mesh; each row of xy is ');
disp(' the coordinates of a point in d-space.              ');

disp('  ');
disp('coordpart(Tmesh,Txy);');
figure(2)
coordpart(Tmesh,Tmeshxy);
disp('cutsize(Tmesh,ans)');
cutsize(Tmesh,ans)
disp(' Hit space to continue ...');
pause;

figure(3)
disp(' ');
disp(' Next is a multilevel method from the "Metis 5.0" package.');
disp(' This will only work if you have Metis and its Matlab interface.');
disp('  ');
disp('metispart(Tmesh,Tmeshxy);');
metispart(Tmesh,Tmeshxy);
disp('cutsize(Tmesh,ans)');
cutsize(Tmesh,ans)
disp('  ');
disp(' Hit space to continue ...');
pause;

figure(3)
disp(' ');
disp(' Next is a hypergrapeh partitioning method from the "hMetis 1.5" package.');
disp(' This will only work if you have hMetis and its Matlab interface.');
disp('  ');
disp('metispart(Tmesh,Tmeshxy);');
hmetispart(Tmesh,Tmeshxy);
disp('cutsize(Tmesh,ans)');
cutsize(Tmesh,ans)
disp('  ');
disp(' Hit space to continue ...');
pause;


figure(4)
disp(' ');
disp(' Next is a multilevel method from the "KaHIP" package.');
disp(' This will only work if you have KaHIP and its Matlab interface.');
disp('  ');
disp('kahippart(Tmesh,Tmeshxy);');
kahippart(Tmesh,Tmeshxy);
disp('cutsize(Tmesh,ans)');
cutsize(Tmesh,ans)
disp('  ');
disp(' Hit space to continue ...');
pause;


figure(5)
disp(' ');
disp(' Next is a multilevel method from the "KaHIP_NE" package.');
disp(' This will only work if you have KaHIP_NE and its Matlab interface.');
disp('  ');
disp('kahipnepart(Tmesh,Tmeshxy);');
kahipnepart(Tmesh,Tmeshxy);
disp('cutsize(Tmesh,ans)');
cutsize(Tmesh,ans)
disp('  ');
disp(' Hit space to continue ...');
pause;



disp(' ');
disp(' Next is the spectral partitioning, which uses the second eigenvector of');
disp(' the Laplacian matrix of the graph, also known as the "Fiedler vector".');
disp('  ');
disp('USIspecpart(Tmesh,Txy);');
figure(6)
USIspecpart(Tmesh,Tmeshxy);
disp('cutsize(Tmesh,ans)');
cutsize(Tmesh,ans)
disp('  ');

disp(' Hit space to continue ...');
pause;

figure(7)
disp(' ');
disp(' Next is inertial partitioning, which uses the coordinates to find');
disp(' a separating line in the plane.');
disp('  ');
disp('inertpart(Tmesh,Tmeshxy);');
USIinertpart(Tmesh,Tmeshxy);

disp('cutsize(Tmesh,ans)');
cutsize(Tmesh,ans)
disp('  ');
disp(' Hit space to continue ...');
pause;
close all;


end % whichdemos(1)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(whichdemos==2) % Vertex separators in graphs.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('          *********************************************')
disp('          ***           Vertex Separators           *** ');
disp('          *********************************************')
disp(' ');

figure(1);
clf reset;
colordef(1,'black')

disp(' ');
disp(' A vertex separator is a set of vertices whose removal partitions the graph.');
disp(' We can convert a vertex partition (which is an edge separator) into a');
disp(' vertex separator by finding a minimum cover with bipartitite matching.');
disp(' ');

disp('part = gspart(Tmesh);');

part = gspart(Tmesh);

disp(' ');
disp(' Here''s the vertex partition (edge separator):');
disp(' ');

disp('gplotpart(Tmesh,Tmeshxy,part)');

gplotpart(Tmesh,Tmeshxy,part)
title('Edge Separator');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' And here''s the vertex separator:');
disp(' ');

disp('sep = vtxsep(Tmesh,part)');

sep = vtxsep(Tmesh,part)

figure(2);
clf reset;
colordef(2,'black')

disp('highlight(Tmesh,Tmeshxy,sep);');

highlight(Tmesh,Tmeshxy,sep);
title('Vertex Separator');
drawnow;

disp(' ');
disp(' How big is the separator?');
disp(' ');

disp('SeparatorSize = length(sep)');

SeparatorSize = length(sep)
xlabel([int2str(SeparatorSize), ' separator vertices'],'visible','on');


disp(' ');
disp(' And how big are the pieces?');
disp(' ');

disp('nonsep = other(sep,Tmesh);');
disp('blocks = components(Tmesh(nonsep,nonsep));');
disp('BlockSizes = hist(blocks,max(blocks))');

nonsep = other(sep,Tmesh);
blocks = components(Tmesh(nonsep,nonsep));
BlockSizes = hist(blocks,max(blocks))

disp(' ');
disp(' Hit space to continue ...');
pause;
disp('  ');

end % whichdemos(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(whichdemos==3) %  Recursive Bisection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('          *********************************************')
disp('          ***          Recursive Bisection          *** ');
disp('          *********************************************')
disp(' ');

figure(1);
clf reset;
colordef(1,'black')

gplotg(Tmesh,Tmeshxy);
nvtx = size(Tmesh,1);
nedge = (nnz(Tmesh)-nvtx)/2;
xlabel([int2str(nvtx) ' vertices, ' int2str(nedge) ' edges'],'visible','on');


disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' Use geometric partitioning to dice the grid into sixteenths.');
disp(' ');

disp('geodice(Tmesh,Tmeshxy);');

geodice(Tmesh,4,Tmeshxy);

disp(' ');
disp(' Hit space to continue ...');
pause;
disp('  ');

end % whichdemos(3)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if any(whichdemos==4) % Nested dissection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' ');
disp('          *********************************************')
disp('          ***           Nested Dissection           *** ');
disp('          *********************************************')
disp(' ');

figure(1);
clf reset;
colordef(1,'black')

disp(' ');
disp(' Nested dissection is a recursively defined sparse matrix permutation');
disp(' based on vertex separators in the graph of a symmetric matrix.');
disp(' Any partitioning algorithm can be used; here we illustrate nested');
disp(' dissection with spectral partitioning.');
disp(' ');
disp(' Here''s the original matrix and its graph...');
disp(' ');

disp('gplotg(Tmesh,Tmeshxy);');

gplotg(Tmesh,Tmeshxy);
title('Finite Element Mesh');
figure(2); 
clf reset;
colordef(2,'black')

disp('spy(Tesh,3);');

spy(Tmesh,3);
title('Coefficient Matrix');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' Here''s the result of factoring the matrix with Gaussian elimination.');

disp(' ');
disp('spy(chol(Smallmesh),3);');

ones(size(Tmesh,1),1);
Tmesh = Tmesh + diag(100*ones(size(Tmesh,1),1));


spy(chol(Tmesh),3);
title('Cholesky Factor');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' Now we use nested dissection to reorder the matrix.');

disp(' ');
disp('nd = specnd(Tmesh);');

nd = specnd(Tmesh);

disp(' ');
disp(' Permute the rows and columns of the matrix with Matlab''s indexing.');

disp(' ');
disp('spy(Tmesh(nd,nd),3);');

figure(1)
spy(Tmesh(nd,nd),3);
title('Coefficient Matrix with Nested Dissection');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' The permuted matrix has a sparser Cholesky factor than the original.');

disp(' ');
disp('spy(chol(Tmesh(nd,nd)),3);');

figure(1);
spy(chol(Tmesh(nd,nd)),3);
title('Cholesky Factor with Nested Dissection');

disp(' ');
disp(' Hit space to continue ...');
pause;

disp(' ');
disp(' The elimination tree, which has one vertex per matrix column,');
disp(' shows the recursive structure of the nested dissection ordering.');
disp(' The height of the tree is related to the complexity of factoring');
disp(' the matrix on a parallel computer.');

disp(' ');
disp('etreeplotg(Tmesh(nd,nd),1);');

figure(2);
etreeplotg(Tmesh(nd,nd),1);
title('Elimination Tree with Nested Dissection');

disp(' ');
disp(' Hit space to continue ...');
pause;

if exist('symamd')  % symamd is new in Matlab 6, use symamd on earlier versions
	disp(' ');
	disp(' Another ordering is minimum degree, which often gives higher trees');
	disp(' but lower fill than nested dissection.  Say "help symamd" for details.');
	disp(' Here''s a summary of the complexity on this matrix:');
	md = symamd(Tmesh);
else
    disp(' ');
	disp(' Another ordering is minimum degree, which often gives higher trees');
	disp(' but lower fill than nested dissection.  Say "help symamd" for details.');
	disp(' Here''s a summary of the complexity on this matrix:');
	md = symamd(Tmesh);
end;

disp(' ');
disp('With no permutation:');
analyze(Tmesh);

disp(' ');
disp('With nested dissection:');
analyze(Tmesh(nd,nd));

disp(' ');
disp('With minimum degree:');
analyze(Tmesh(md,md));

disp(' ');
disp(' Hit space to continue ...');
pause;

end % whichdemos(4)



end % loop over meshes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final notes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(' ');
disp(' ');
disp(' ');
disp('  **********************************************************************')
disp('  *** For more information, see meshpart.html or say "help meshpart" ***');
disp('  **********************************************************************')
disp(' ');

format;
