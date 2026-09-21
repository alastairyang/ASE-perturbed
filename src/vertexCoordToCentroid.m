function coord_interp = vertexCoordToCentroid(coord, index)
    % Get the coordinates defined on the element centroid
    % 
    % example: x_elem = VertexCoordToCentroid(md.mesh.x, md.mesh.elements2d);
    %.         y_elem = VertexCoordToCentroid(md.mesh.y, md.mesh.elements2d);
    summation = [1;1;1];
    coord_interp = coord(index)*summation/3; % data is defined on vertices; take the average
end
