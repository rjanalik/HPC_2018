
function writeMeshToVTKFile(prefix, ElementList, PointList, ElementParams, PointParams, type)

numberOfPoints   = size(PointList, 1);
numberOfElements = size(ElementList, 1);
numberOfVertices = size(ElementList, 2);

% 2. read the .node file
% %%%%%%%%%%%%%%%%%%%%%%

% this opens a file in text 't mode for read access 'r'
filename = strcat(prefix, '.vtk');

fprintf('writting mesh file %s\n', filename);

fout = fopen(filename, 'w');
fprintf(fout,'# vtk DataFile Version 5.10\n');
fprintf(fout,'Hexahedral mesh with data\n');
fprintf(fout,'ASCII\n');
fprintf(fout,'DATASET UNSTRUCTURED_GRID\n');
fprintf(fout,'POINTS %d float\n', numberOfPoints);


% now write the PointList
% -----------------------


for i = 1:numberOfPoints
    x_i = PointList(i,:);
    fprintf(fout,'%25.16e %25.16e %25.16e\n', x_i(1), x_i(2), x_i(3));
    
end

fprintf(fout,'\n');
entries = (numberOfVertices+1)*numberOfElements;
fprintf(fout,'CELLS %d %d\n', numberOfElements, entries);
first_number = min(ElementList(:));

for e = 1:numberOfElements
    
    v_e = ElementList(e, :);
    
    v_e = v_e - first_number;
    
    fprintf(fout,'%d ', numberOfVertices);
    
    for i=1:numberOfVertices
        fprintf(fout,'%d ', v_e(i));
    end
    
    fprintf(fout, '\n');
    
end


fprintf(fout,'\n');
fprintf(fout,'CELL_TYPES %d\n', numberOfElements);

for e = 1:numberOfElements
    
    fprintf(fout,'%d\n', type);
    
end

fprintf(fout,'\n');

nsets = size(ElementParams,2);
if (nsets > 0)
    
    fprintf(fout,'CELL_DATA %d\n', numberOfElements);
    fprintf(fout,'SCALARS Surface float 1\n');
    fprintf(fout,'LOOKUP_TABLE default\n');
    for n = 1:numberOfElements
        
        fprintf(fout,'%25.16e\n', ElementParams(n,1));
        
    end
end



nsets = size(PointParams,2);
if (nsets > 0)
    
    fprintf(fout,'POINT_DATA %d\n', numberOfPoints);
    fprintf(fout,'SCALARS NodePartitioning float 1\n');
    fprintf(fout,'LOOKUP_TABLE default\n');
    for n = 1:numberOfPoints
        
        fprintf(fout,'%25.16e\n', PointParams(n,1));
        
    end
end


fclose(fout);
end


