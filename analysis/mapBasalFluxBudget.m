% compute the volumetric strain heating
n_glen    = 3;
secinyear = 24*3600*365;

% --- Strain heating ---
B       = md.materials.rheology_B;
A       = B.^(-n_glen);
eps_dot = md.results.TransientSolution(end).StrainRateeffective;   % s^-1, 90280x1
Q       = 2 * A.^(-1/n_glen) .* eps_dot.^((n_glen+1)/n_glen);     % W/m^3, 90280x1

% --- Geometry ---
z = md.mesh.z;                                                      % 90280x1 [m]
H = md.results.TransientSolution(end).Thickness;                   % 90280x1 [m]

% --- Vertical velocity ---
% Vz in ISSM is positive upward -> flip sign so positive = downward
mda   = augmentModel(md);
w_avg = mda.plotDepthAverage(md.results.TransientSolution(end).Vz, 0);  % 18056x1 [m/yr]
w_avg = -w_avg / secinyear;                                         % m/s, positive downward

% --- Material properties ---
k   = md.materials.thermalconductivity;
rho = md.materials.rho_ice;
cp  = md.materials.heatcapacity;

% --- Compute ---
dT_base = strain_heating_basal_temp_3d(z, Q, H, w_avg, k, rho, cp);

fprintf('dT_base range: %.4f to %.4f K\n', min(dT_base), max(dT_base));

%%
i = 16859;
zr_col = zrel_3d(:, i);
q_col  = Q_3d(:, i);
Hi     = H_col(i);
ai     = alpha(i);

fprintf('ai = %.4e, Hi = %.2f, ai*Hi = %.2f\n', ai, Hi, ai*Hi);

exponent = -ai * (Hi - zr_col);
fprintf('exponents: '); fprintf('%.2f ', exponent); fprintf('\n');

exp_vals = exp(exponent);
fprintf('exp vals: '); fprintf('%.4e ', exp_vals); fprintf('\n');

G = (1 - exp_vals) / (ai * k);
fprintf('G vals: '); fprintf('%.4e ', G); fprintf('\n');

dT_i = trapz(zr_col, q_col .* G);
fprintf('dT for node %d: %.6f K\n', i, dT_i);



%%

dT_base = strain_heating_basal_temp_3d(z, Q, H, w_avg, k, rho, cp);

fprintf('dT_base range: %.4f to %.4f K\n', min(dT_base), max(dT_base));
fprintf('dT_base mean:  %.4f K\n',          mean(dT_base));

% Visual check
figure;
scatter(md.mesh.x2d, md.mesh.y2d, 10, dT_base, 'filled');
clim([0,20])
colorbar; colormap(hot);
title('\Delta T_{base} from strain heating [K]');
axis equal tight;




%%
function dT_base = strain_heating_basal_temp_3d(z, Q, H, w_avg, k, rho, cp)

    N = numel(w_avg);
    M = numel(z) / N;

    % Reshape to M x N (row 1 = bed, row M = surface)
    z_3d = reshape(z, N, M)';    % M x N
    Q_3d = reshape(Q, N, M)';    % M x N

    % Bed elevation = layer 1 (row 1)
    z_bed = z_3d(1, :);          % 1 x N  [m absolute elevation]

    % Height above bed for each node
    zrel_3d = z_3d - z_bed;      % M x N  [m], 0 at bed, H at surface

    % Ice thickness = top - bottom elevation
    H_col = zrel_3d(M, :);       % 1 x N  [m], should match md.results...Thickness

    alpha = (rho * cp / k) * w_avg(:)';   % 1 x N  [m^-1]

    dT_base = zeros(N, 1);

    for i = 1:N
        zr_col = zrel_3d(:, i);   % M x 1, height above bed
        q_col  = Q_3d(:, i);      % M x 1
        Hi     = H_col(i);        % scalar, ice thickness
        ai     = alpha(i);

        % Numerically stable Green's function
        % G(0,z') = (1 - exp(-alpha*(H-z'))) / (alpha*k)
        %         = (H - z') / k   for alpha -> 0
        if abs(ai) < 1e-10
            G = (Hi - zr_col) / k;
        else
            G = (1 - exp(-ai * (Hi - zr_col))) / (ai * k);
        end

        dT_base(i) = trapz(zr_col, q_col .* G);
    end
end




function dz = compute_layer_weights(z_mat, N, M)
% Trapezoidal weights for non-uniform vertical layers
% z_mat: N x M, monotonically increasing along M dimension
% Returns dz: N x M weights for integration

    dz = zeros(N, M);

    % Interior layers: central difference weights
    dz(:, 2:M-1) = 0.5 * (z_mat(:, 3:M) - z_mat(:, 1:M-2));

    % End layers: one-sided
    dz(:, 1) = 0.5 * (z_mat(:, 2) - z_mat(:, 1));
    dz(:, M) = 0.5 * (z_mat(:, M) - z_mat(:, M-1));

end
