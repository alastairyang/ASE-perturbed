% function runmeFuncAISThermal(name, sim_steps, Ts_spec, smb_spec, sliding_spec, array_id, nproc)
% %RUNMEFUNCTHERMAL thermal modeling of AIS
% %   Ts_spec : "max" or "min". Max or min surface temperature from several 
% %            regional climate model output
% %   smb_spec: "max" or "min". Max or min surface mass balance from several
% %            regional climate model outputs
% %   sliding_spec: "softsliding" or "hardsliding" changes exponent m in
% %            Budd's sliding parameterization

% ----------------------------------------
% -------------- PARAMETERS --------------
sim_steps = [7,8];
name = 'ASE';
Ts_spec = "max";
smb_spec = "max";
sliding_spec = "softsliding";
nproc   = 20; 
t_final = 275; % duration of perturbation run. (through 2300)
t_phase1 = 1;  % duration of constant time stepping (relaxing). 


% -----------------------------------------

% ---------- name - expfile pair
drainage_filename = 'gh_thwaites_pi_boundary_subsampled_20perc.exp';

% ---------- END OF PARAMETERS --------------
% get the environmental variables
common_data_dir  = getenv("COMMON_DATA");
project_data_dir = getenv("PROJECT_DATA");
result_data_dir  = getenv("PROJECT_RESULTS");

if isempty(common_data_dir)
    error("!!! The file directories are empty -- run the startup script first!")
end

% set project name
modelername = 'Yang';
model_name = [name '_perturbed'];

% dataset paths
% make folder for experiment (i.e., this combination of parameters)
foldername = [result_data_dir 'model'];
if ~exist(foldername, 'dir')
    mkdir(foldername)
end

% model name prefix
model_identifier  = string(modelername)+"_"+string(Ts_spec)+"_"+string(smb_spec)+"_"+string(sliding_spec);
model_name_prefix = string(name)+ "_"+model_identifier+"_";
% convert to char
model_identifier  = convertStringsToChars(model_identifier);
model_name_prefix = convertStringsToChars(model_name_prefix);

% file paths
switch sliding_spec
    case "hardsliding"
        par_path = [project_data_dir 'par-files/thwaites-PIG-hardsliding.par'];
    case "softsliding"
        par_path = [project_data_dir 'par-files/thwaites-PIG-softsliding.par'];
    otherwise
        error('Unknow sliding specification!')
end
vel_path         = [common_data_dir 'velocity/antarctica_ice_velocity_450m_v2.nc'];
GHF_path         = [project_data_dir 'saved-data/GHF-forcing/sequential-gaussian-sim/scaled/'];
eliza_model_path = [common_data_dir 'other-simulations/Dawson_Antarctica_Initialization.mat'];
save_folder      = [result_data_dir 'thermal-training-data/Thwaites-PIG/training/raw/'];


md.cluster = generic('name',oshostname(),'np',nproc);
md.cluster

disp(['Detected hostname: ', oshostname()])

% step iteration
for steps = sim_steps

    % Run steps
    org=organizer('repository',foldername,'prefix',model_name_prefix,'steps',steps);

    if perform(org, 'Meshing') % {{{1  STEP 1 % meshing the model domain
        % load the drainage boundary exp file

        md = model();
        md = triangle(md, [project_data_dir 'exp-files/' drainage_filename],2000);
        md.miscellaneous.name = model_name;

        % set up cluster information
        md.settings.waitonlock = 0;
        md.cluster = generic('name',oshostname(),'np',nproc);

        % load ice velocity
        vx = double(ncread(vel_path, 'VX'));
        vy = double(ncread(vel_path, 'VY'));
        x  = double(ncread(vel_path, 'x'));
        y  = double(ncread(vel_path, 'y'));

        vx_md = InterpFromGridToMesh(x, flipud(y), flipud(vx'), ...
            md.mesh.x, md.mesh.y, 0);
        vy_md = InterpFromGridToMesh(x, flipud(y), flipud(vy'), ...
            md.mesh.x, md.mesh.y, 0);

        % material prop for getting strain rate
        md.materials.rheology_B = cuffey(265.15)*ones(md.mesh.numberofvertices,1); % -23 degree depth averaged
        md.materials.rheology_n = 3*ones(md.mesh.numberofelements,1);
        md.materials.rheology_law = 'CuffeyTemperate';

        % for thermal model: use strain rate as a proxy
        % to refine the mesh
        md = mechanicalproperties(md,vx_md,vy_md);
        eps_eff = md.results.strainrate.effectivevalue;
        % from element to vertices
        elem_x = vertexCoordToCentroid(md.mesh.x, md.mesh.elements);
        elem_y = vertexCoordToCentroid(md.mesh.y, md.mesh.elements);
        F_eps_eff = scatteredInterpolant(elem_x, elem_y, eps_eff);
        eps_eff_vertice = F_eps_eff(md.mesh.x, md.mesh.y);

        md = bamg(md, 'field', 3e3*eps_eff_vertice, ...
                     'hmax', 20000, 'hmin', 4000,'err',2); % original: hmax=40000, err = 3; to be replaced with hmax=20000, err=0.5
    
        clear vel_md vx_md vy_md 
            
        % specify the number of processor based of 1:1000 element ratio
        np = min(round(md.mesh.numberofelements/1000), feature('numcores'));
        cluster = generic('name', oshostname(), 'np', np);
        cluster.interactive = 1;

        % solver and tolerance
        md.stressbalance.abstol=NaN;
        md.stressbalance.restol=1e-4;
        md.stressbalance.maxiter=50;
        % md.toolkits.DefaultAnalysis = bcgslbjacobioptions(); % just for the first 2D inversion 
        md.settings.solver_residue_threshold = 1e-3; % I set it; Denis didn't.

        savemodel(org,md);
    end

    if perform(org, 'Parameterization') % {{{2 step 2 Parameterize model

        md = loadmodel(org, 'Meshing');
        md = parameterize(md, par_path);

        % md = extrude(md,8,2.5);
        % md = setflowequation(md,'HO','all'); 
        md = setflowequation(md,'SSA','all'); 

        % set up cluster information
        md.settings.waitonlock = 0;
        md.cluster = generic('name',oshostname(),'np',nproc);

        savemodel(org,md)

    end

    if perform(org, 'Inversion') % {{{3 step 3 inversion
        stepname = 'Inversion';
        md = loadmodel(org, 'Parameterization');

        % set up cluster information
        md.settings.waitonlock = 0;
        md.cluster = generic('name',oshostname(),'np',nproc);

        % inversion velocity
        md.inversion.vx_obs = md.initialization.vx;
        md.inversion.vy_obs = md.initialization.vy;
        md.inversion.vz_obs = md.initialization.vz;
        md.inversion.vel_obs = md.initialization.vel;
        
        md.inversion = m1qn3inversion(md.inversion);
        md.inversion.iscontrol=1;
        md.verbose = verbose('solution',false,'control',true);
        md.transient.amr_frequency = 0;
        
        % Cost functions
        md.inversion.cost_functions = [101 103 501];
        md.inversion.cost_functions_coefficients = zeros(md.mesh.numberofvertices,numel(md.inversion.cost_functions));
        
        % fast flow
        pos = md.inversion.vel_obs > 200;
        switch sliding_spec
            case "hardsliding"
                md.inversion.maxiter = 50;
                md.inversion.cost_functions_coefficients(:,1) = 1000; % 800
                md.inversion.cost_functions_coefficients(:,2) = 0.1; % 1
                md.inversion.cost_functions_coefficients(:,3) = 5e-6; % 1e-6
                md.inversion.min_parameters = 0.01*ones(md.mesh.numberofvertices,1);
                md.inversion.max_parameters = 1e5*ones(md.mesh.numberofvertices,1);
                pos=find(md.mask.ice_levelset>0);
                md.inversion.cost_functions_coefficients(pos,1)=0;
            case "softsliding"
                md.inversion.maxiter = 100;
                md.inversion.cost_functions_coefficients(:,1) = 50;
                % md.inversion.cost_functions_coefficients(pos,1) = 500; % 100
                md.inversion.cost_functions_coefficients(:,2) = 0.1; % 1
                md.inversion.cost_functions_coefficients(:,3) = 2.0448e-5; % 2.0448e-5;
                md.inversion.min_parameters = 0.01*ones(md.mesh.numberofvertices,1);
                md.inversion.max_parameters = 1e8*ones(md.mesh.numberofvertices,1);
                pos=find(md.mask.ice_levelset>0);
                md.inversion.cost_functions_coefficients(pos,1)=0;            
        end
        % pos=find(md.mask.ice_levelset>0);
        % md.inversion.cost_functions_coefficients(pos,1)=0;
        
        % Controls
        md.inversion.control_parameters = {'FrictionCoefficient'};
        md.inversion.maxsteps = 100;
        md.inversion.control_scaling_factors = 1;
        md.inversion.dxmin = 0.01; % 0.01
        md.inversion.gttol = 1e-8;
        % Additional parameters
        md.stressbalance.restol = 0.001;
        md.stressbalance.reltol = 0.1;
        md.stressbalance.abstol = NaN;
        % 
        % % you can make the following nan for inversion if it's really
        % % struggling (e.g. high exponent in the sliding
        % % parameterization), but be sure to turn this back on 
        % md.settings.solver_residue_threshold = nan;
        % 
        disp('Starting inversion!') 
        md = solve(md,'sb');
        md = loadresultsfromcluster(md);

        if isa(md.mesh, 'mesh2d') % haven't extruded into 3d model
            disp("   Extruding into 3D model...")
            md = extrude(md, 5, 2.5);
            md = setflowequation(md, 'HO','all');
            
            md.friction.coefficient = repmat(md.results.StressbalanceSolution.FrictionCoefficient,...
                                             md.mesh.numberoflayers,1);
            md.initialization.vx    = repmat(md.results.StressbalanceSolution.Vx,...
                                             md.mesh.numberoflayers,1);
            md.initialization.vy    = repmat(md.results.StressbalanceSolution.Vy,...
                                             md.mesh.numberoflayers,1);
            md.initialization.vel   = sqrt(md.initialization.vx.^2 + ...
                                         md.initialization.vy.^2);
            md.initialization.pressure = repmat(md.results.StressbalanceSolution.Pressure,...
                                                md.mesh.numberoflayers,1);

            % do a stress balance solve with the 3D model to get vz
            md.inversion.iscontrol=0;
            md = solve(md,'sb');
            md = loadresultsfromcluster(md);

            md.initialization.vx = md.results.StressbalanceSolution.Vx;
            md.initialization.vy = md.results.StressbalanceSolution.Vy;
            md.initialization.vz = md.results.StressbalanceSolution.Vz;
            md.initialization.vel = sqrt(md.initialization.vx.^2 + ...
                                         md.initialization.vy.^2 + ...
                                         md.initialization.vz.^2);
            md.initialization.pressure=md.results.StressbalanceSolution.Pressure;

        else
            md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient;
            md.initialization.vx = md.results.StressbalanceSolution.Vx;
            md.initialization.vy = md.results.StressbalanceSolution.Vy;
            md.initialization.vz = md.results.StressbalanceSolution.Vz;
            md.initialization.vel = sqrt(md.initialization.vx.^2 + ...
                                         md.initialization.vy.^2 + ...
                                         md.initialization.vz.^2);
            md.initialization.pressure=md.results.StressbalanceSolution.Pressure;
        end
        % turn back on
        md.settings.solver_residue_threshold = 1e-4;

        savemodel(org,md);
    end % }}}

    if perform(org, 'Inversion3D') % {{{4 step 4 inversion with 3D model
        stepname = 'Inversion3D';

        md = loadmodel(org, 'Inversion');

        md.inversion.iscontrol = 1;
        md.inversion.maxiter = 10;
        md = solve(md,'sb');
        md = loadresultsfromcluster(md);

        % update the fields
        md.friction.coefficient = md.results.StressbalanceSolution.FrictionCoefficient;
        md.initialization.vx = md.results.StressbalanceSolution.Vx;
        md.initialization.vy = md.results.StressbalanceSolution.Vy;
        md.initialization.vz = md.results.StressbalanceSolution.Vz;
        md.initialization.vel = sqrt(md.initialization.vx.^2 + ...
                                     md.initialization.vy.^2 + ...
                                     md.initialization.vz.^2);
        md.initialization.pressure=md.results.StressbalanceSolution.Pressure;

        savemodel(org,md);
    end
        

    if perform(org, 'PrescribeThermalBC') % {{{5 step 5 prescribe thermal boundary
        stepname = 'PrescribeThermalBC';

        md = loadmodel(org, 'Inversion3D');
        md.inversion.iscontrol = 0;

        % set up cluster information
        md.settings.waitonlock = 0;
        md.cluster = generic('name',oshostname(),'np',nproc);

        % load continental model for temperature profile at influx
        % boundary
        md_eliza = load(eliza_model_path).md;
        md2d_e = model();
        md2d_e.mesh.x = md_eliza.mesh.x2d;
        md2d_e.mesh.y = md_eliza.mesh.y2d;
        md2d_e.mesh.elements = md_eliza.mesh.elements2d;
        md2d_e.mesh.numberofelements = md_eliza.mesh.numberofelements2d;
        md2d_e.mesh.numberofvertices = md_eliza.mesh.numberofvertices2d;
        disp('----Eliza model is loaded!')

        % my model but 2D
        md2d = model();

        md2d.mesh.x = md.mesh.x2d;
        md2d.mesh.y = md.mesh.y2d;
        md2d.mesh.elements = md.mesh.elements2d;
        md2d.mesh.numberofelements = md.mesh.numberofelements2d;
        md2d.mesh.numberofvertices = md.mesh.numberofvertices2d;

        wf_val = 1e-6; % initial water fraction
        wf_2d  = wf_val*ones(md.mesh.numberofvertices2d,1); 

        switch Ts_spec
            case "max"
                load([project_data_dir 'Ts-RCM/Ts_max.mat']);
                func_ts = scatteredInterpolant(Ts_max_struct.x, ...
                                           Ts_max_struct.y,...
                                           Ts_max_struct.Ts);
                % md.miscellaneous.name = [md.miscellaneous.name '; Ts: max'];
            case "min"
                load([project_data_dir 'Ts-RCM/Ts_min.mat']);
                func_ts = scatteredInterpolant(Ts_min_struct.x, ...
                                           Ts_min_struct.y,...
                                           Ts_min_struct.Ts);
                % md.miscellaneous.name = [md.miscellaneous.name '; Ts: min'];
            otherwise
                error("Unknown spec")
        end
        switch smb_spec
            case "max"
                load([project_data_dir 'SMB-RCM/smb_max.mat']);
                func_smb = scatteredInterpolant(smb_max_struct.x, ...
                                           smb_max_struct.y,...
                                           smb_max_struct.Ts);
                % md.miscellaneous.name = [md.miscellaneous.name '; SMB: max'];
            case "min"
                load([project_data_dir 'SMB-RCM/smb_min.mat']);
                func_smb = scatteredInterpolant(smb_min_struct.x, ...
                                           smb_min_struct.y,...
                                           smb_min_struct.Ts);
                % md.miscellaneous.name = [md.miscellaneous.name '; SMB: min'];
            otherwise
                error("Unknown spec")
        end
        smb = func_smb(md2d.mesh.x, md2d.mesh.y);
        Ts  = func_ts(md2d.mesh.x,  md2d.mesh.y);

        % get ice surface temperature and water fraction from Eliza's model
        Tinterp  = griddata(md_eliza.mesh.x, md_eliza.mesh.y, ...
                            md_eliza.mesh.z, ...
                            md_eliza.results.ThermalSolution.Temperature, ...
                            md.mesh.x, md.mesh.y,md.mesh.z,'nearest');
        
        STmean = -40 + 273.15;
        BTmean = -10 + 273.15; % use only as a first guess for model initialization
        
        % thermal model: initial and surface boundary conditions
        md.thermal.spctemperature = NaN*ones(md.mesh.numberofvertices, 1);
        md.thermal.spctemperature(find(md.mesh.vertexonsurface)) = Ts;
        md.thermal.spctemperature(find(md.mesh.vertexonboundary)) = Tinterp(find(md.mesh.vertexonboundary));
        md.initialization.temperature = NaN*ones(md.mesh.numberofvertices, 1);
        md.initialization.temperature(find(md.mesh.vertexonsurface)) = Ts;
        % make sure that no point exceeds pressure melting point
        replicate=repmat(md.geometry.surface-md.mesh.z,1,size(md.thermal.spctemperature,2));
        apm = md.materials.meltingpoint-md.materials.beta*md.materials.rho_ice*md.constants.g*replicate+1e-5;
        md.thermal.spctemperature(find(md.thermal.spctemperature > apm)) = apm(md.thermal.spctemperature > apm);
        
        % initialize temperature
        % linear interpolate the englacial temperature
        dz = md.mesh.z(find(md.mesh.vertexonbase)) - md.mesh.z(find(md.mesh.vertexonsurface));
        dTdz = (BTmean - STmean)./dz;
        dTdz_HO = repmat(dTdz, md.mesh.numberoflayers,1);
        BTmean_HO = repmat(BTmean, md.mesh.numberoflayers*numel(dz), 1);
        md.initialization.temperature = BTmean_HO + dTdz_HO.*(md.mesh.z - md.geometry.base);
        md.initialization.temperature(isnan(md.initialization.temperature)) = STmean;

        % specify surface mass balance
        md.smb.mass_balance(find(md.mesh.vertexonsurface)) = smb;
        
        % initialize water (only required if using enthalpy method)
        md.initialization.watercolumn = zeros(md.mesh.numberofvertices,1);
        md.initialization.waterfraction = zeros(md.mesh.numberofvertices,1);
        md.initialization.waterfraction(find(md.mesh.vertexonbase)) = wf_2d;
        md.initialization.waterfraction(md.initialization.waterfraction<=0) = 1e-6; % non-zero minima

        % save the model
        savemodel(org,md);
    end

    if perform(org,'ThermalSteadyState') % {{{6 step 6 thermal steady state model
        % here we do initialize the rheology and do a transient
        % relaxation
        stepname = 'ThermalSteadyState';
        md = loadmodel(org, 'PrescribeThermalBC');

        % ensure we are using default analysis for HO model
        md2 = model();
        md.toolkits.DefaultAnalysis = md2.toolkits.DefaultAnalysis;
        clear md2

        % set up cluster information
        md.settings.waitonlock = 0;

        md.verbose.solution = 1;
        md.verbose.convergence = 0;
        md.verbose.control = 0;

        md.inversion.iscontrol = 0;
        md = solve(md, 'sb');
        md = loadresultsfromcluster(md);

        % update velocity field
        md.initialization.vx = md.results.StressbalanceSolution.Vx;
        md.initialization.vy = md.results.StressbalanceSolution.Vy;
        md.initialization.vz = md.results.StressbalanceSolution.Vz;
        md.initialization.vel = sqrt(md.initialization.vx.^2 + ...
            md.initialization.vy.^2 + ...
            md.initialization.vz.^2);

        % geothermal heat flux
        % load GHF posterior mean from the inference
        posterior_path = '/data-archive/ASE-inference/';
        ghf_posterior = load([posterior_path 'GHF/posterior_low_mean_high.mat']);
        ghf_mean_interp = InterpFromGridToMesh(ghf_posterior.x, ghf_posterior.y,...
                                               ghf_posterior.posterior_mean,...
                                               md.mesh.x2d, md.mesh.y2d, nan);
        % slightly extrapolate to areas under ice shelf
        % GHF there doesn't matter anyways
        nonnan_idx = ~isnan(ghf_mean_interp);
        nonnan_x   = md.mesh.x2d(nonnan_idx);
        nonnan_y   = md.mesh.y2d(nonnan_idx);
        nonnan_ghf = ghf_mean_interp(nonnan_idx);
        ghf_mean_interp_func = scatteredInterpolant(nonnan_x, nonnan_y,...
                                                    nonnan_ghf, 'linear','linear');
        ghf_mean_interp_full = ghf_mean_interp_func(md.mesh.x2d, md.mesh.y2d) * 1e-3; % convert to W/m^2
        
        ghf_prefactor = 0.7; % ad hoc measure to make interior cooler

        md.basalforcings.geothermalflux = ghf_prefactor * repmat(ghf_mean_interp_full, md.mesh.numberoflayers,1);

        % solver specification
        md.timestepping.time_step = 0;
        md.timestepping.final_time = 0;
        md.thermal.isenthalpy = 1;
        md.thermal.stabilization = 3;
        md.thermal.isdynamicbasalspc = 1;
        md.thermal.fe = 'P1xP2';
        md.hydrology.spcwatercolumn = NaN;
        md.thermal.maxiter = 200;
        md.verbose.solution = 1;
        md.verbose.convergence = 1;
        md.verbose.control = 0;

        md = solve(md,'thermal');
        md = loadresultsfromcluster(md);

        % update ice rheology
        md.materials.rheology_B = cuffey(md.results.ThermalSolution.Temperature);

        % ----------- do a long-term relaxation ------------------------
        % first a fixed timestep to ensure stability
        md.transient.isthermal = 0;
        md.transient.isgroundingline = 1;

        md.timestepping.time_step  = 0.005;
        md.timestepping.final_time = 3;
        
        % necessary once isgroundingline is turned on 
        md.calving.calvingrate = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.meltingrate = zeros(md.mesh.numberofvertices,1);
        md.frontalforcings.ablationrate = zeros(md.mesh.numberofvertices,1);

        md.settings.output_frequency = 200; 
        md = solve(md,'tr');
        md = loadresultsfromcluster(md);

        % save an immediate version just in case
        save([result_data_dir 'model/intermediate/ASE_before_thermal_transient_relax.mat'],'md')

        % --------- relaxed with adaptive time stepping ------------------
        % coupled with thermal solution
        md = transientrestart(md);

        % with the relaxed geometry, update the boundary conditon
        replicate=repmat(md.geometry.surface-md.mesh.z,1,size(md.thermal.spctemperature,2));
        apm = md.materials.meltingpoint-md.materials.beta*md.materials.rho_ice*md.constants.g*replicate+1e-5;
        md.thermal.spctemperature(find(md.thermal.spctemperature > apm)) = apm(md.thermal.spctemperature > apm);
        md.transient.isthermal = 1;

        md.materials.rheology_law = 'Cuffey';
        md.initialization.temperature   = md.results.ThermalSolution.Temperature;
        md.initialization.waterfraction = md.results.ThermalSolution.Waterfraction;
        md.initialization.watercolumn   = md.results.ThermalSolution.Watercolumn;
        
        md.timestepping = timesteppingadaptive();
        md.timestepping.time_step_min = 0.005;
        md.timestepping.time_step_max = 0.5;
        md.timestepping.final_time = md.timestepping.start_time + 100;

        md.settings.output_frequency = 100; 
        
        % solve
        md = solve(md,'tr');
        md = loadresultsfromcluster(md);

        savemodel(org,md);
    end

    if perform(org,'ShelfCollapseConstantStep') % {{{7 step 7 perturbation: shelf collapse - Small constant time stepping
        % here we completely remove the ice shelf
        stepname = 'ShelfCollapseConstantStep';
        md = loadmodel(org, 'ThermalSteadyState');
        md = transientrestart(md);

        % whole sim duration
        t_phase2 = t_final;

        % with the relaxed geometry, update the thermal dirichlet boundary conditon
        replicate=repmat(md.geometry.surface-md.mesh.z,1,size(md.thermal.spctemperature,2));
        apm = md.materials.meltingpoint-md.materials.beta*md.materials.rho_ice*md.constants.g*replicate+1e-5;
        md.thermal.spctemperature(find(md.thermal.spctemperature > apm)) = apm(md.thermal.spctemperature > apm);

        % regular time stepping to avoid blowing up
        % run for small number of years
        md.timestepping = timestepping();
        md.timestepping.time_step = 0.01;
        md.timestepping.final_time = t_phase1;
        md.settings.output_frequency = 20;


        % update the mask
        shelf_idx = md.mask.ocean_levelset(1:md.mesh.numberofvertices2d) < 0;
        init_ice_mask = md.mask.ice_levelset(1:md.mesh.numberofvertices2d);
        collapse_ice = init_ice_mask;
        collapse_ice(shelf_idx) = 1;

        n_layer = md.mesh.numberoflayers;

        init_ice_mask = repmat(init_ice_mask, n_layer, 1);
        collapse_ice_mask = repmat(collapse_ice, n_layer, 1); % collapse after one timestep
        phase1_ice_mask    = repmat(collapse_ice, n_layer,1); % hold position
        phase2_ice_mask    = repmat(collapse_ice, n_layer,1); % hold position

        % add time stamp
        init_ice_mask_t = [init_ice_mask; md.timestepping.start_time];
        collapse_ice_mask_t = [collapse_ice_mask; md.timestepping.time_step];
        phase1_ice_mask_t   = [phase1_ice_mask; t_phase1];
        phase2_ice_mask_t   = [phase2_ice_mask; t_phase2];

        md.levelset.spclevelset = [init_ice_mask_t, collapse_ice_mask_t, phase1_ice_mask_t, phase2_ice_mask_t];
        
        % transient settings
        md.transient.ismovingfront = 1;

        % solve
        md = solve(md,'tr');
        md = loadresultsfromcluster(md);        

        savemodel(org, md);
    end

    if perform(org, 'ShelfCollapseAdaptiveRegularSliding') % {{{8 step 8 perturbation: shelf collapse - adaptive time stepping 
        md = loadmodel(org, 'ShelfCollapseConstantStep');
        md = transientrestart(md);

        md.settings.output_frequency = 50;
        md.stressbalance.requested_outputs = {'default', 'StrainRateeffective'}; 
        
        md.timestepping = timesteppingadaptive();
        md.timestepping.start_time = t_phase1;
        md.timestepping.final_time = md.timestepping.start_time + t_final;

        % solve
        md = solve(md,'tr');
        md = loadresultsfromcluster(md);        


        savemodel(org, md);


    end

    if perform(org, 'ShelfCollapseAdaptiveTempSliding') % {{{9 step 9 perturbation: shelf collapse - adaptive time stepping and temperature-dependent sliding coefficient
        md = loadmodel(org, 'ShelfCollapseConstantStep');
        md = transientrestart(md);

        % create an empty struct field to collect 
        results = [];

        md.settings.output_frequency = 4;
        md.stressbalance.requested_outputs = {'default', 'StrainRateeffective'}; 
        
        md.timestepping = timesteppingadaptive();

        delta_t = 1; % update rheology every 1 year
        nt = floor(t_final / delta_t);
        
        for ii = 1:nt
            disp(['Starting run segment no. ' num2str(ii) ' !'])
            if ii ~= 1
                md.materials.rheology_B = cuffey(md.results.TransientSolution(end).Temperature);
                
                results = [results, md.results.TransientSolution];

                end_time_prev = md.timestepping.final_time;
                md = transientrestart(md);
            else
                end_time_prev = t_phase1;
            end

            % update temperature B.C.
            replicate=repmat(md.geometry.surface-md.mesh.z,1,size(md.thermal.spctemperature,2));
            apm = md.materials.meltingpoint-md.materials.beta*md.materials.rho_ice*md.constants.g*replicate+1e-5;
            md.thermal.spctemperature(find(md.thermal.spctemperature > apm)) = apm(md.thermal.spctemperature > apm);


            md.timestepping.start_time = end_time_prev;
            md.timestepping.final_time = md.timestepping.start_time + delta_t;

            % solve
            md = solve(md,'tr');
            md = loadresultsfromcluster(md);        
        end

        results = [results, md.results.TransientSolution];

        % subsample to 1/10
        results = results(1:10:end); % ~ 60 rows 
        save("perturbation_results_all.mat", "results")

        savemodel(org, md);


    end

        
end

