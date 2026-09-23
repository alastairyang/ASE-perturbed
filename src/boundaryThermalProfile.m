function T_func = boundaryThermalProfile(md, ghf, Ts, smb, solver)
%BOUNDARYTHERMALPROFILE create temperature profiles at the influx boundary 

    % get the boundary vertices at the base
    basebound_idx = md.mesh.vertexonbase & md.mesh.vertexonboundary;
    surfbound_idx = md.mesh.vertexonsurface & md.mesh.vertexonboundary;
    xb = md.mesh.x(basebound_idx);
    yb = md.mesh.y(basebound_idx);
    zs = md.mesh.z(surfbound_idx);
    zb = md.mesh.z(basebound_idx);

    H_bc = zs - zb;

    if numel(ghf) > md.mesh.numberofvertices2d; error("Expect 2d GHF input"); end
    if numel(Ts)  > md.mesh.numberofvertices2d; error("Expect 2d Ts input"); end
    if numel(smb) > md.mesh.numberofvertices2d; error("Expect 2d SMB input"); end

    GHF_bc = ghf;
    smb_bc = smb;
    Ts_bc  = Ts;
    
    % get temperature at influx boundaries
    % Parameters you do not need to change
    dsdx = 0.0001; % surface slope (only useful for non-ice-shelf; turned off)
    u = 0; % m/a, along-flow velocity
    bm = 0; 
    dz_ref = 10; % m, spacing
    dTdx = 0; % horizontal temperature gradient

    % create a thermal mode
    z_vec = cell(numel(H_bc),1);
    T_vec = cell(numel(H_bc),1);
    for jj = 1:numel(H_bc)
        if H_bc(jj) <= dz_ref % if thickness is smaller than ref interval
            dz = H_bc(jj)*0.2;
        else
            dz = dz_ref;
        end
        tmd = thermal1d(H_bc(jj), dsdx, dz, smb_bc(jj), bm, Ts_bc(jj), GHF_bc(jj), u, dTdx, "DJ", solver);
        tmd = tmd.solveNonlinearGI();

        z_vec{jj} = tmd.sol.x;
        T_vec{jj} = tmd.sol.y(1,:);
    end

    % construct coords
    coords_md = [];
    temp_md = [];
    for jj = 1:numel(H_bc)
        % collect coordinates
        xv = repmat(xb(jj), numel(z_vec{jj}),1);
        yv = repmat(yb(jj), numel(z_vec{jj}),1);
        coords_temp = [xv, yv, zb(jj) + transpose(z_vec{jj})];
        coords_md = [coords_md; coords_temp];

        % collect temperature
        temp_md = [temp_md; transpose(T_vec{jj})];
    end

    % concat all coordinates and temperature to create an
    % interpolant
    T_func = scatteredInterpolant(coords_md(:,1),...
                                         coords_md(:,2),...
                                         coords_md(:,3),...
                                         temp_md);
    
end