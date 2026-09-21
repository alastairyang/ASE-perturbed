function [px, py, pmag] = drivingstressFromResults(md, time_i)
%DRIVINGSTRESS_FROM_RESULTS calculate the driving stress, but not from
%md.geometry; from md.results.TransientSolution.
%
%   NOTE: The outputs are defined on the element centroid. Use
%   vertexCoordToCentroid to get the element centroid coordinate; you can
%   then use scatteredInterpolant for interpolating onto any grid / mesh. 
%
%   Input:
%       md: ISSM model class
%       time_i: time index / No. row in solution table.

    % taken from ISSM built-in function "slope.m"
    if dimension(md.mesh)==2
	    numberofelements=md.mesh.numberofelements;
	    numberofnodes=md.mesh.numberofvertices;
	    index=md.mesh.elements;
	    x=md.mesh.x; y=md.mesh.y;

    	surf=md.results.TransientSolution(time_i).Surface;
        thickness = md.results.TransientSolution(time_i).Thickness;

    else % if 3D model
	    numberofelements=md.mesh.numberofelements2d;
	    numberofnodes=md.mesh.numberofvertices2d;
	    index=md.mesh.elements2d;
	    x=md.mesh.x2d; y=md.mesh.y2d;

        top_layer = md.mesh.numberoflayers;

        surf = project2d(md, md.results.TransientSolution(time_i).Surface, top_layer);
        thickness = project2d(md, md.results.TransientSolution(time_i).Thickness, top_layer);
    end

    %compute nodal functions coefficients N(x,y)=alpha x + beta y + gamma
    [alpha, beta]=GetNodalFunctionsCoeff(index,x,y);
    
    summation=[1;1;1];
    sx=(surf(index).*alpha)*summation;
    sy=(surf(index).*beta)*summation;
    s=sqrt(sx.^2+sy.^2);

    % Average thickness over elements
    thickness_bar=(thickness(index(:,1))+thickness(index(:,2))+thickness(index(:,3)))/3;
    
    % get driving stress for x and y components
    px=-md.materials.rho_ice*md.constants.g*thickness_bar.*sx;
    py=-md.materials.rho_ice*md.constants.g*thickness_bar.*sy;
    pmag=sqrt(px.^2+py.^2);
    
    
end

