%% Loading
md_path = 'results/model/ASE_Yang_max_max_softsliding_ShelfCollapseAdaptiveRegularSliding.mat';
load(md_path);

% % load results
% load("perturbation_results_all.mat") % "results"
% md.results.TransientSolution = results; 
%%

thwaites_fl = expread('data/exp-files/thwaites-flowline.exp');
pig_fl      = expread('data/exp-files/PIG-flowline.exp');
thwaites_fl_xy = unique([thwaites_fl.x, thwaites_fl.y], "rows");
pig_fl_xy      = unique([pig_fl.x,      pig_fl.y], "rows");

thwaites_fl = interpCurviLine(thwaites_fl_xy, 500);
pig_fl      = interpCurviLine(pig_fl_xy     , 500);

colormap_path = '/home/donglaiyang/Documents/Georgia-Tech/Research/common-data-set/ScientificColourMaps8/';
load([colormap_path 'bukavu/bukavu.mat'])
load([colormap_path 'vik/vik.mat'])
%% extract ice thickness profile, driving stress profile, and basal temperature
n_time = size(md.results.TransientSolution, 2);

taud_thwaites = cell(n_time, 1);
taud_pig      = cell(n_time, 1);
Tb_thwaites   = cell(n_time, 1);
Tb_pig        = cell(n_time, 1);
pmp_thwaites  = cell(n_time, 1);
pmp_pig       = cell(n_time, 1);
vel_thwaites  = cell(n_time, 1);
vel_pig       = cell(n_time, 1);

x_element = vertexCoordToCentroid(md.mesh.x2d, md.mesh.elements2d);
y_element = vertexCoordToCentroid(md.mesh.y2d, md.mesh.elements2d);

for ii = 1:n_time
    % driving stress (element-centroid interpolant)
    [~, ~, taud] = drivingstressFromResults(md, ii);
    taud_func = scatteredInterpolant(x_element, y_element, taud);
    taud_thwaites{ii} = taud_func(thwaites_fl(:,1), thwaites_fl(:,2));
    taud_pig{ii}      = taud_func(pig_fl(:,1),      pig_fl(:,2));

    % basal temperature (vertex interpolant, layer 1 = base)
    Tb = project2d(md, md.results.TransientSolution(ii).Temperature, 1);
    Tb_func = scatteredInterpolant(md.mesh.x2d, md.mesh.y2d, Tb);
    Tb_thwaites{ii} = Tb_func(thwaites_fl(:,1), thwaites_fl(:,2));
    Tb_pig{ii}      = Tb_func(pig_fl(:,1),      pig_fl(:,2));

    % pressure melting point
    h   = project2d(md, md.results.TransientSolution(ii).Thickness, 1);
    pmp = md.materials.meltingpoint - md.materials.beta * ...
          md.materials.rho_ice * md.constants.g * h + 1e-5;
    pmp_func = scatteredInterpolant(md.mesh.x2d, md.mesh.y2d, pmp);
    pmp_thwaites{ii} = pmp_func(thwaites_fl(:,1), thwaites_fl(:,2));
    pmp_pig{ii}      = pmp_func(pig_fl(:,1),      pig_fl(:,2));

    % velocity magnitude (vertex interpolant, layer 1 = base)
    vel = project2d(md, md.results.TransientSolution(ii).Vel, 1);
    vel_func = scatteredInterpolant(md.mesh.x2d, md.mesh.y2d, vel);
    vel_thwaites{ii} = vel_func(thwaites_fl(:,1), thwaites_fl(:,2));
    vel_pig{ii}      = vel_func(pig_fl(:,1),      pig_fl(:,2));
end

% bed topography (static — outside loop)
bed = project2d(md, md.geometry.bed, 1);
bed_func = scatteredInterpolant(md.mesh.x2d, md.mesh.y2d, bed);
bed_thwaites = bed_func(thwaites_fl(:,1), thwaites_fl(:,2));
bed_pig      = bed_func(pig_fl(:,1),      pig_fl(:,2));


%% Make a 2D map showing where the flowlines are
mda = augmentModel(md);
md2d = mda.create_2D_mesh();
% get baseline velocity vector
vx = project2d(md, md.initialization.vx, md.mesh.numberoflayers);
vy = project2d(md, md.initialization.vy, md.mesh.numberoflayers);

% ── Subsample nodes for quiver ────────────────────────────────────────────
quiver_skip = 10;
idx_q = 1:quiver_skip:length(md.mesh.x2d);

xq  = md.mesh.x2d(idx_q);
yq  = md.mesh.y2d(idx_q);
vxq = vx(idx_q);
vyq = vy(idx_q);

% ── Log-scale the arrow lengths ───────────────────────────────────────────
speed     = sqrt(vxq.^2 + vyq.^2);
log_speed = log10(speed + 1);
unit_vx   = vxq ./ (speed + eps);
unit_vy   = vyq ./ (speed + eps);

arrow_scale = 8000;
vx_scaled   = unit_vx .* log_speed * arrow_scale;
vy_scaled   = unit_vy .* log_speed * arrow_scale;

% ── Plot ──────────────────────────────────────────────────────────────────
figure('Position',[1000,500,1000,1000])
patch('Faces',    md.mesh.elements2d, ...
      'Vertices', [md.mesh.x2d md.mesh.y2d], ...
      'FaceVertexCData', bed, ...
      'FaceColor', 'interp', 'EdgeColor', 'none');
colormap(bukavu);
clim([-2000,2000])

% ── Colorbar: wider + bigger font ─────────────────────────────────────────
cb = colorbar;
cb.FontSize = 24;                  % 3× a typical 12pt default
cb.Position(2) = cb.Position(2) + 0.1;
cb.Position(3) = cb.Position(3) * 2;   % triple the width
cb.Position(4) = cb.Position(4) * 0.8;
cb.Label.String = 'Bed topography (m)';
hold on;

quiver(xq, yq, vx_scaled, vy_scaled, 0, ...
       'Color', [0.2 0.2 0.2], ...
       'LineWidth', 0.8, ...
       'MaxHeadSize', 0.3);

plot(thwaites_fl(:,1), thwaites_fl(:,2), '-', 'LineWidth', 4, 'Color', 'blue');
plot(pig_fl(:,1),      pig_fl(:,2),      '-', 'LineWidth', 4, 'Color', 'red');

scatter(thwaites_fl(1,1), thwaites_fl(1,2), 300, 'b', 'filled');
scatter(pig_fl(1,1),      pig_fl(1,2),      300, 'r', 'filled');

hold off;
axis off;

exportgraphics(gcf,'figs/bedmap_flowlines.png','Resolution',300)

%% make tiledlayout plot
plot_time_interval = 4;
fs = 26;

figure('Position', [1000, 500, 1400, 1100]);
tiledlayout(3, 2, 'TileSpacing', 'tight', 'Padding', 'compact');

% ── Precompute distances ──────────────────────────────────────────────────
dx = diff(thwaites_fl(:,1)); dy = diff(thwaites_fl(:,2));
dist_thwaites = cumsum([0; sqrt(dx.^2 + dy.^2)]);   

dx = diff(pig_fl(:,1));      dy = diff(pig_fl(:,2));
dist_pig = cumsum([0; sqrt(dx.^2 + dy.^2)]);

epsilon = 0.5;

% ═════════════════════════════════════════════════════════════════════════
% ROW 1 — Driving Stress
% ═════════════════════════════════════════════════════════════════════════

% ── [1,1] Thwaites – Driving Stress ──────────────────────────────────────
nexttile(1);
tau_base = taud_thwaites{1};

yyaxis left
plot(dist_thwaites, tau_base, 'k-', 'LineWidth', 2);
ylabel('\tau_d (Pa)  {\color{black}\bf—}');
ax = gca; ax.YAxis(1).Color = 'k';

yyaxis right
hold on;
for i = 2:plot_time_interval:n_time
    delta_tau = taud_thwaites{i} - tau_base;
    pos = delta_tau; pos(pos < 0) = 0;
    fill([dist_thwaites; flipud(dist_thwaites)], [pos; zeros(size(pos))], ...
         'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    neg = delta_tau; neg(neg > 0) = 0;
    fill([dist_thwaites; flipud(dist_thwaites)], [neg; zeros(size(neg))], ...
         'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
hold off;
ylabel('\Delta\tau_d (Pa)');
ax.YAxis(2).Color = 'k';
title('Thwaites'); xticklabels([]);
grid on; box on; set(gca, 'FontSize', fs);

% ── [1,2] PIG – Driving Stress ───────────────────────────────────────────
nexttile(2);
tau_base = taud_pig{1};

yyaxis left
plot(dist_pig, tau_base, 'k-', 'LineWidth', 2);
ylabel('\tau_d (Pa)  {\color{black}\bf—}');
ax = gca; ax.YAxis(1).Color = 'k';

yyaxis right
hold on;
for i = 2:plot_time_interval:n_time
    delta_tau = taud_pig{i} - tau_base;
    pos = delta_tau; pos(pos < 0) = 0;
    fill([dist_pig; flipud(dist_pig)], [pos; zeros(size(pos))], ...
         'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    neg = delta_tau; neg(neg > 0) = 0;
    fill([dist_pig; flipud(dist_pig)], [neg; zeros(size(neg))], ...
         'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
hold off;
ylabel('\Delta\tau_d (Pa)');
ax.YAxis(2).Color = 'k';
title('PIG'); xticklabels([]);
grid on; box on; set(gca, 'FontSize', fs);

% ═════════════════════════════════════════════════════════════════════════
% ROW 2 — Velocity
% ═════════════════════════════════════════════════════════════════════════

% ── [2,1] Thwaites – Velocity ─────────────────────────────────────────────
nexttile(3);
vel_base = vel_thwaites{1};

yyaxis left
plot(dist_thwaites, vel_base, 'k-', 'LineWidth', 2);
ylabel('Velocity (m/yr)  {\color{black}\bf—}');
ax = gca; ax.YAxis(1).Color = 'k';

yyaxis right
hold on;
for i = 2:plot_time_interval:n_time
    delta_vel = vel_thwaites{i} - vel_base;
    pos = delta_vel; pos(pos < 0) = 0;
    fill([dist_thwaites; flipud(dist_thwaites)], [pos; zeros(size(pos))], ...
         'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    neg = delta_vel; neg(neg > 0) = 0;
    fill([dist_thwaites; flipud(dist_thwaites)], [neg; zeros(size(neg))], ...
         'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
hold off;
ylabel('\DeltaV (m/yr)');
ax.YAxis(2).Color = 'k';
xticklabels([]);
grid on; box on; set(gca, 'FontSize', fs);

% ── [2,2] PIG – Velocity ──────────────────────────────────────────────────
nexttile(4);
vel_base = vel_pig{1};

yyaxis left
plot(dist_pig, vel_base, 'k-', 'LineWidth', 2);
ylabel('Velocity (m/yr)  {\color{black}\bf—}');
ax = gca; ax.YAxis(1).Color = 'k';

yyaxis right
hold on;
for i = 2:plot_time_interval:n_time
    delta_vel = vel_pig{i} - vel_base;
    pos = delta_vel; pos(pos < 0) = 0;
    fill([dist_pig; flipud(dist_pig)], [pos; zeros(size(pos))], ...
         'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    neg = delta_vel; neg(neg > 0) = 0;
    fill([dist_pig; flipud(dist_pig)], [neg; zeros(size(neg))], ...
         'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
hold off;
ylabel('\DeltaV (m/yr)');
ax.YAxis(2).Color = 'k';
xticklabels([]);
grid on; box on; set(gca, 'FontSize', fs);

% ═════════════════════════════════════════════════════════════════════════
% ROW 3 — Bed & ΔT_b
% ═════════════════════════════════════════════════════════════════════════

% ── [3,1] Thwaites – Bed & ΔT_b ──────────────────────────────────────────
nexttile(5);
Tb_base        = Tb_thwaites{1};
pmp_base       = pmp_thwaites{1};
already_thawed = abs(Tb_base - pmp_base) <= epsilon;

yyaxis left
plot(dist_thwaites, bed_thwaites, 'k-', 'LineWidth', 2);
ylabel('Bed elev. (m)  {\color{black}\bf—}');
ax = gca; ax.YAxis(1).Color = 'k';

yyaxis right
hold on;
for i = 2:plot_time_interval:n_time
    delta_Tb = Tb_thwaites{i} - Tb_base;
    pos = delta_Tb; pos(pos < 0) = 0;
    fill([dist_thwaites; flipud(dist_thwaites)], [pos; zeros(size(pos))], ...
         'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    neg = delta_Tb; neg(neg > 0) = 0;
    fill([dist_thwaites; flipud(dist_thwaites)], [neg; zeros(size(neg))], ...
         'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
yl = ylim;
thawed_y = nan(size(dist_thwaites));
thawed_y(already_thawed) = yl(2);
plot(dist_thwaites, thawed_y, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 12);
hold off;
ylabel('\DeltaT_b (K)  [{\color[rgb]{0.85,0.33,0.10}\bf— thawed}]');
ax.YAxis(2).Color = 'k';
xlabel('Distance from grounding line (m)');
grid on; box on; set(gca, 'FontSize', fs);

% ── [3,2] PIG – Bed & ΔT_b ───────────────────────────────────────────────
nexttile(6);
Tb_base        = Tb_pig{1};
pmp_base       = pmp_pig{1};
already_thawed = abs(Tb_base - pmp_base) <= epsilon;

yyaxis left
plot(dist_pig, bed_pig, 'k-', 'LineWidth', 2);
ylabel('Bed elev. (m)  {\color{black}\bf—}');
ax = gca; ax.YAxis(1).Color = 'k';

yyaxis right
hold on;
for i = 2:plot_time_interval:n_time
    delta_Tb = Tb_pig{i} - Tb_base;
    pos = delta_Tb; pos(pos < 0) = 0;
    fill([dist_pig; flipud(dist_pig)], [pos; zeros(size(pos))], ...
         'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    neg = delta_Tb; neg(neg > 0) = 0;
    fill([dist_pig; flipud(dist_pig)], [neg; zeros(size(neg))], ...
         'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none');
end
yl = ylim;
thawed_y = nan(size(dist_pig));
thawed_y(already_thawed) = yl(2);
plot(dist_pig, thawed_y, '-', 'Color', [0.85 0.33 0.10], 'LineWidth', 12);
hold off;
ylabel('\DeltaT_b (K)  [{\color[rgb]{0.85,0.33,0.10}\bf— thawed}]');
ax.YAxis(2).Color = 'k';
xlabel('Distance from grounding line (m)');
grid on; box on; set(gca, 'FontSize', fs);

drawnow;

%% plot animated cross section
% ── Settings ─────────────────────────────────────────────────────────────
% gif_filename = 'thwaites_temp_xsection.gif';
frame_delay  = 0.15;
fs           = 20;
resolution   = [1000, 50];
n_time       = size(md.results.TransientSolution, 2);
thwaites_exp = 'data/exp-files/thwaites-flowline.exp';
pig_exp      = 'data/exp-files/PIG-flowline.exp';

exp_choice = pig_exp;
fl_choice  = pig_fl;
gif_filename = 'pig_temp_xsection.gif';


% ── Regular grid for continuous plotting ──────────────────────────────────
n_s    = 300;
n_zeta = 100;
zeta_grid = linspace(0, 1, n_zeta);

% ── Precompute flowline distance ──────────────────────────────────────────
dx = diff(pig_fl(:,1));
dy = diff(pig_fl(:,2));
dist_fl_km = cumsum([0; sqrt(dx.^2 + dy.^2)]) / 1e3;
n_fl = length(dist_fl_km);

% ── Basal shear stress (static) ───────────────────────────────────────────
[~, ~, bs] = basalstress(md);
bs_func = scatteredInterpolant(md.mesh.x, md.mesh.y, bs, 'linear', 'nearest');
bs_fl   = bs_func(pig_fl(:,1), pig_fl(:,2));   % Pa

% ── Pre-scan ΔT range ─────────────────────────────────────────────────────
fprintf('Pre-scanning ΔT range...\n');
[~,~,~,~,~, T_0] = SectionValues(md, ...
    md.results.TransientSolution(1).Temperature, exp_choice, resolution);

dT_all = [];
for ii = 2:n_time
    [~,~,~,~,~, T_curr] = SectionValues(md, ...
        md.results.TransientSolution(ii).Temperature, exp_choice, resolution);
    dT_all = [dT_all; T_curr(:) - T_0(:)];
end
clim_val = max(abs(prctile(dT_all, [2 98])));
clim_val = max(clim_val, 0.01);
fprintf('ΔT color range: ±%.3f K\n', clim_val);

% ── Pre-scan frictional heating ───────────────────────────────────────────
fprintf('Pre-scanning frictional heating range...\n');
fric_all = zeros(n_fl, n_time-1);
for ii = 2:n_time
    vel_base = project2d(md, md.results.TransientSolution(ii).Vel, 1);
    vel_func = scatteredInterpolant(md.mesh.x2d, md.mesh.y2d, vel_base, 'linear', 'nearest');
    vel_fl   = vel_func(pig_fl(:,1), pig_fl(:,2));   % m/yr
    fric_all(:, ii-1) = bs_fl .* (vel_fl / md.constants.yts);  % W/m²
end

fric_0        = fric_all(:, 1);                          % reference = first transient step
dfric_all     = fric_all - fric_0;                       % change from initial
fric_ylim     = [0, prctile(fric_all(:), 99)];
dfric_clim    = max(abs(prctile(dfric_all(:), [2 98])));
dfric_clim    = max(dfric_clim, 1e-4);

% ── Line colors (time progression) ────────────────────────────────────────
line_colors = parula(n_time-1);

% ── Build GIF ─────────────────────────────────────────────────────────────
fig = figure('Position', [1000, 200, 1600, 500], 'Color', 'w');
t_init = md.results.TransientSolution(1).time;
t_end  = md.results.TransientSolution(n_time).time;

for ii = 2:n_time
    sol = md.results.TransientSolution(ii);

    % ── Temperature section ───────────────────────────────────────────────
    [~,~,~, z_curr, s_curr, T_curr] = SectionValues(md, sol.Temperature, exp_choice, resolution);
    [~,~,~,~,~,      S_curr]        = SectionValues(md, sol.Surface,     exp_choice, resolution);
    [~,~,~,~,~,      B_curr]        = SectionValues(md, sol.Base,        exp_choice, resolution);

    thickness = S_curr - B_curr;
    zeta_pts  = (z_curr - B_curr) ./ max(thickness, 1);
    zeta_pts  = max(0, min(1, zeta_pts));
    s_pts     = s_curr / 1e3;
    dT        = T_curr - T_0;

    s_grid_vec        = linspace(min(s_pts), max(s_pts), n_s);
    [S_grid, Z_grid]  = meshgrid(s_grid_vec, zeta_grid);
    T_reg  = fliplr(griddata(s_pts, zeta_pts, T_curr, S_grid, Z_grid, 'linear'));
    dT_reg = fliplr(griddata(s_pts, zeta_pts, dT,     S_grid, Z_grid, 'linear'));

    fric_curr  = fric_all(:,  ii-1);
    dfric_curr = dfric_all(:, ii-1);
    t_curr     = sol.time;

    clf;

    % ── [1,1] Absolute temperature ────────────────────────────────────────
    subplot(2, 2, 1);
    pcolor(S_grid, Z_grid, T_reg);
    shading interp;
    colormap(gca, parula(256));
    cb1 = colorbar;
    cb1.Label.String   = 'Temperature (K)';
    cb1.Label.FontSize = fs - 2;
    ylabel('\zeta  (0=base, 1=surface)', 'FontSize', fs);
    title(sprintf('Temperature   t = %.1f yr', t_curr), 'FontSize', fs);
    set(gca, 'FontSize', fs);
    xticklabels([]);
    xlim([min(s_pts) max(s_pts)]); ylim([0 1]);

    % ── [1,2] ΔT from initial ─────────────────────────────────────────────
    subplot(2, 2, 2);
    pcolor(S_grid, Z_grid, dT_reg);
    shading interp;
    colormap(gca, vik);
    clim([-clim_val, clim_val]);
    cb2 = colorbar;
    cb2.Label.String   = sprintf('\\DeltaT (K) from t = %.1f yr', t_init);
    cb2.Label.FontSize = fs - 2;
    ylabel('\zeta  (0=base, 1=surface)', 'FontSize', fs);
    title(sprintf('\\DeltaT   step %d / %d', ii-1, n_time-1), 'FontSize', fs);
    set(gca, 'FontSize', fs);
    xticklabels([]);
    xlim([min(s_pts) max(s_pts)]); ylim([0 1]);

    % ── [2,1] Absolute frictional heating — accumulating lines ────────────
    subplot(2, 2, 3);
    hold on;
    for jj = 1:ii-2
        plot(dist_fl_km, fric_all(:, jj), ...
             'Color', [line_colors(jj,:), 0.4], 'LineWidth', 1.0);
    end
    plot(dist_fl_km, fric_curr, ...
         'Color', line_colors(ii-1,:), 'LineWidth', 2.5);
    hold off;
    ylim(fric_ylim);
    xlim([0 max(dist_fl_km)]);
    xlabel('Distance from grounding line (km)', 'FontSize', fs);
    ylabel('W m^{-2}', 'FontSize', fs);
    title('\tau_b \cdot u_b  (absolute)', 'FontSize', fs);
    set(gca, 'FontSize', fs);
    box on; grid on;
    cb3 = colorbar;
    colormap(gca, parula(n_time-1));
    clim([t_init, t_end]);
    cb3.Label.String   = 'Time (yr)';
    cb3.Label.FontSize = fs - 2;

    % ── [2,2] Δ frictional heating from initial — accumulating lines ──────
    subplot(2, 2, 4);
    hold on;
    for jj = 1:ii-2
        plot(dist_fl_km, dfric_all(:, jj) * 1e3, ...
             'Color', [line_colors(jj,:), 0.4], 'LineWidth', 1.0);
    end
    plot(dist_fl_km, dfric_curr * 1e3, ...
         'Color', line_colors(ii-1,:), 'LineWidth', 2.5);
    yline(0, 'k--', 'LineWidth', 1.0);
    hold off;
    ylim([0, 50]);
    xlim([0 max(dist_fl_km)]);
    xlabel('Distance from grounding line (km)', 'FontSize', fs);
    ylabel('mW m^{-2}', 'FontSize', fs);
    title('\Delta(\tau_b \cdot u_b)  from t = ' + string(t_init) + ' yr', 'FontSize', fs);
    set(gca, 'FontSize', fs);
    box on; grid on;
    cb4 = colorbar;
    colormap(gca, parula(n_time-1));
    clim([t_init, t_end]);
    cb4.Label.String   = 'Time (yr)';
    cb4.Label.FontSize = fs - 2;

    drawnow;

    % ── Write frame ───────────────────────────────────────────────────────
    frame = getframe(fig);
    [imind, cm] = rgb2ind(frame2im(frame), 256);
    if ii == 2
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', Inf, 'DelayTime', frame_delay);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
    end

    fprintf('Frame %d / %d\n', ii-1, n_time-1);
end
fprintf('Done → %s\n', gif_filename);


%% Animated 2D map — bed topo + velocity change quiver + ice front + grounding line
gif_filename = 'figs/thwaites_map_animation.gif';
frame_delay  = 0.15;
fs           = 18;
n_time       = size(md.results.TransientSolution, 2);

% ── Static background: bed topography ────────────────────────────────────
bed = project2d(md, md.geometry.bed, 1);

% ── Subsampled quiver indices (static) ────────────────────────────────────
quiver_skip = 10;
idx_q = 1:quiver_skip:length(md.mesh.x2d);
xq = md.mesh.x2d(idx_q);
yq = md.mesh.y2d(idx_q);

% ── Reference velocity at t=1 ─────────────────────────────────────────────
sol_0  = md.results.TransientSolution(1);
vx_0   = project2d(md, sol_0.Vx, md.mesh.numberoflayers);
vy_0   = project2d(md, sol_0.Vy, md.mesh.numberoflayers);
vx_0q  = vx_0(idx_q);
vy_0q  = vy_0(idx_q);

% ── Pre-scan Δvelocity range for consistent arrow scaling ─────────────────
fprintf('Pre-scanning delta-velocity range...\n');
max_dspeed_all = 0;
for ii = 2:n_time
    sol    = md.results.TransientSolution(ii);
    vx_ii  = project2d(md, sol.Vx, md.mesh.numberoflayers);
    vy_ii  = project2d(md, sol.Vy, md.mesh.numberoflayers);
    dvx    = vx_ii(idx_q) - vx_0q;
    dvy    = vy_ii(idx_q) - vy_0q;
    dspd   = sqrt(dvx.^2 + dvy.^2);
    max_dspeed_all = max(max_dspeed_all, prctile(dspd, 99));  % robust max
end
arrow_scale = 8000;
fprintf('Max |Δv| (99th pct): %.1f m/yr\n', max_dspeed_all);

% ── Precompute isoline edges once ─────────────────────────────────────────
fprintf('Precomputing isoline edges...\n');
mda  = augmentModel(md);
md2d = mda.create_2D_mesh();
[~, edges] = isoline(md2d.md2d, ...
    project2d(md, md.results.TransientSolution(1).MaskOceanLevelset, 1), ...
    'value', 0);

% ── Build GIF (start from ii=2) ───────────────────────────────────────────
fig = figure('Position', [1200, 200, 1000, 1000], 'Color', 'w');

for ii = 2:n_time
    sol    = md.results.TransientSolution(ii);
    t_curr = sol.time;
    t_init = sol_0.time;

    % ── Δvelocity from first timestep ─────────────────────────────────────
    vx_ii  = project2d(md, sol.Vx, md.mesh.numberoflayers);
    vy_ii  = project2d(md, sol.Vy, md.mesh.numberoflayers);
    dvx    = vx_ii(idx_q) - vx_0q;
    dvy    = vy_ii(idx_q) - vy_0q;
    dspeed = sqrt(dvx.^2 + dvy.^2);

    % log-scale arrow lengths, preserve direction of Δv
    log_dspeed = log10(dspeed + 1);
    vx_scaled  = (dvx ./ (dspeed + eps)) .* log_dspeed * arrow_scale;
    vy_scaled  = (dvy ./ (dspeed + eps)) .* log_dspeed * arrow_scale;

    % ── Isolines ──────────────────────────────────────────────────────────
    gl = isoline(md2d.md2d, project2d(md, sol.MaskOceanLevelset, 1), ...
                 'value', 0, 'output', 'matrix', 'edges', edges);
    icefront = isoline(md2d.md2d, project2d(md, sol.MaskIceLevelset, 1), ...
                       'value', 0, 'output', 'matrix', 'edges', edges);

    clf;

    % ── Bed topography patch ──────────────────────────────────────────────
    patch('Faces',           md.mesh.elements2d, ...
          'Vertices',        [md.mesh.x2d, md.mesh.y2d], ...
          'FaceVertexCData', bed, ...
          'FaceColor',       'interp', ...
          'EdgeColor',       'none', ...
          'DisplayName',     'Bed topography');
    colormap(gca, bukavu);
    clim([-2000, 2000]);

    cb = colorbar;
    cb.FontSize       = fs;
    cb.Position(2)    = cb.Position(2) + 0.10;
    cb.Position(3)    = cb.Position(3) * 2;
    cb.Position(4)    = cb.Position(4) * 0.80;
    cb.Label.String   = 'Bed topography (m)';
    cb.Label.FontSize = fs;

    hold on;

    % ── Δvelocity quiver ──────────────────────────────────────────────────
    quiver(xq, yq, vx_scaled, vy_scaled, 0, ...
           'Color',       [0.15 0.15 0.15], ...
           'LineWidth',   0.8, ...
           'MaxHeadSize', 0.3, ...
           'DisplayName', sprintf('\\Deltav from t=%.1f yr', t_init));

    % ── Grounding line ────────────────────────────────────────────────────
    plot(gl(:,1), gl(:,2), '-', ...
         'Color',       [0.9 0.5 0.0], ...
         'LineWidth',   2.5, ...
         'DisplayName', 'Grounding line');

    % ── Ice front ─────────────────────────────────────────────────────────
    plot(icefront(:,1), icefront(:,2), '-', ...
         'Color',       'black', ...
         'LineWidth',   2.5, ...
         'DisplayName', 'Ice front');

    hold off;

    % ── Legend + title ────────────────────────────────────────────────────
    legend('Location', 'northwest', 'FontSize', fs-2, 'TextColor', 'w', ...
           'Color', [0.15 0.15 0.15]);
    title(sprintf('\\Deltav from t = %.1f yr  →  t = %.1f yr', t_init, t_curr), ...
          'FontSize', fs+2, 'FontWeight', 'bold');
    axis off;

    drawnow;

    % ── Write frame ───────────────────────────────────────────────────────
    frame = getframe(fig);
    [imind, cm] = rgb2ind(frame2im(frame), 256);
    if ii == 2
        imwrite(imind, cm, gif_filename, 'gif', 'Loopcount', Inf, 'DelayTime', frame_delay);
    else
        imwrite(imind, cm, gif_filename, 'gif', 'WriteMode', 'append', 'DelayTime', frame_delay);
    end

    fprintf('Frame %d / %d   (t = %.1f yr)\n', ii-1, n_time-1, t_curr);
end
fprintf('Done → %s\n', gif_filename);


