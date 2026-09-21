%% Basal-temperature sensitivity to thickness and vertical velocity
% z = 0: bed; z = H: surface.
% w > 0: upward motion; w < 0: downward motion.
%
% Governing equation:
%   k*Tzz - rho*cp*w*Tz + qdot = 0
%
% Boundary conditions:
%   T(H) = Ts
%   -k*Tz(0) = gb
%
% Main output:
%   dTb(i,j) = basal warming from strain heating alone [K]
%              for thickness H(i) and velocity wYr(j).
%
% Heating is prescribed, not coupled to temperature or ice deformation.

clear; close all; clc;

%% Physical parameters: SI units internally
rho = 917;                  % Ice density [kg/m^3]
cp  = 2009;                 % Specific heat capacity [J/(kg K)]
k   = 2.1;                  % Thermal conductivity [W/(m K)]
kappa = k/(rho*cp);          % Thermal diffusivity [m^2/s]

secondsPerYear = 365.25*24*3600;

% Prescribed basal-localized strain heating
QA  = 0.03;                 % Integrated heating [W/m^2]
ell = 100;                  % Heating decay scale above bed [m]

% Background boundary conditions: do not affect dTb
Ts = -25;                   % Surface temperature [deg C]
gb = 0.05;                  % Upward basal conductive flux [W/m^2]

%% Parameter sweep
H   = linspace(300,3000,91);       % Ice thickness [m]
wYr = linspace(-0.05,0.05,101);    % Vertical velocity [m/yr]
w   = wYr/secondsPerYear;         % Vertical velocity [m/s]

nH = numel(H);
nW = numel(w);
Nz = 2001;                       % Vertical quadrature points

dTb = zeros(nH,nW);
Tbg = zeros(nH,nW);

%% Evaluate Green-function integrals
for i = 1:nH
    Hi = H(i);
    z = linspace(0,Hi,Nz).';      % Column vector, height above bed

    % Exponential heating profile, normalized to integral QA
    normFactor = -expm1(-Hi/ell); % Stable evaluation of 1-exp(-Hi/ell)
    qdot = QA/(ell*normFactor) * exp(-z/ell);  % [W/m^3]

    % Local Peclet number: Nz-by-nW array
    Pe = ((Hi-z)/kappa) * w;

    % Basal Green kernel [m^2 K/W]
    % Implicit expansion multiplies each column by (Hi-z)/k.
    Kb = ((Hi-z)/k) .* expRatio(Pe);

    % Heating-induced basal warming [K]
    dTb(i,:) = trapz(z, qdot .* Kb, 1);

    % Background basal temperature without strain heating [deg C]
    PeH = Hi*w/kappa;
    Tbg(i,:) = Ts + gb*(Hi/k)*expRatio(PeH);
end

% Formal total basal temperature for the prescribed linear problem
Tb = Tbg + dTb;

%% Check the zero-velocity solution against an analytical result
[~,jZero] = min(abs(wYr));

% For the normalized exponential heating profile:
% dTb(H,0) = QA/k * [H - ell + H/(exp(H/ell)-1)]
dTbDiffusion = (QA/k) * ...
    (H(:) - ell + H(:)./expm1(H(:)/ell));

relativeError = max(abs(dTb(:,jZero)-dTbDiffusion) ...
                    ./dTbDiffusion);

fprintf('Thermal diffusivity: %.2f m^2/yr\n', ...
        kappa*secondsPerYear);
fprintf('Zero-velocity quadrature check: %.2e relative error\n', ...
        relativeError);
fprintf('Integrated strain heating: %.1f mW/m^2\n',1000*QA);

%% Figure: heating profile and sensitivities
figure('Color','w','Position',[80 80 1200 850]);

% ---------------------------------------------------------------
% Panel 1: basal-localized volumetric heating
% ---------------------------------------------------------------
subplot(2,2,1);
hold on;

Hshow = [300 1000 3000];

for h = Hshow
    zp = linspace(0,h,Nz).';
    qp = QA/(ell*(-expm1(-h/ell))) * exp(-zp/ell);

    plot(1000*qp,zp,'LineWidth',1.8, ...
        'DisplayName',sprintf('H = %g m',h));
end

ylim([0 500]);
xlabel('Volumetric strain heating [mW m^{-3}]');
ylabel('Height above bed [m]');
title('Heating concentrated near the bed');
legend('Location','northeast');
grid on; box on;

% ---------------------------------------------------------------
% Panel 2: thickness sensitivity at selected velocities
% ---------------------------------------------------------------
subplot(2,2,2);
hold on;

wShow = [-0.05 -0.025 0 0.025 0.05];

for v = wShow
    [~,j] = min(abs(wYr-v));

    semilogy(H/1000,dTb(:,j),'LineWidth',1.8, ...
        'DisplayName',sprintf('w = %+.3f m/yr',wYr(j)));
end

% Explicit setting also ensures logarithmic axes with hold on
set(gca,'YScale','log');
xlabel('Ice thickness [km]');
ylabel('Strain-heating basal warming [K]');
title('Thickness sensitivity');
legend('Location','northwest');
grid on; box on;

% ---------------------------------------------------------------
% Panel 3: velocity sensitivity at selected thicknesses
% ---------------------------------------------------------------
subplot(2,2,3);
hold on;

Hshow = [600 1500 3000];

for h = Hshow
    [~,i] = min(abs(H-h));

    semilogy(wYr,dTb(i,:),'LineWidth',1.8, ...
        'DisplayName',sprintf('H = %.0f m',H(i)));
end

set(gca,'YScale','log');
xlabel('Vertical velocity [m/yr]: downward (-), upward (+)');
ylabel('Strain-heating basal warming [K]');
title('Velocity sensitivity through zero');
legend('Location','northwest');
grid on; box on;

% ---------------------------------------------------------------
% Panel 4: combined thickness-velocity sensitivity
% ---------------------------------------------------------------
subplot(2,2,4);

contourf(H/1000,wYr,log10(dTb.'),25,'LineColor','none');
hold on;
plot([H(1) H(end)]/1000,[0 0],'k--','LineWidth',1.2);

xlabel('Ice thickness [km]');
ylabel('Vertical velocity [m/yr]');
title('Combined sensitivity');
cb = colorbar;
cb.Label.String = 'log_{10}(basal warming / 1 K)';
box on;

%% Example numerical output at H = 1500 m
[~,iExample] = min(abs(H-1500));
jExample = [1 jZero nW];

fprintf('\nExample at H = %.0f m:\n',H(iExample));

exampleTable = table( ...
    wYr(jExample).', ...
    dTb(iExample,jExample).', ...
    Tbg(iExample,jExample).', ...
    Tb(iExample,jExample).', ...
    'VariableNames', ...
    {'w_m_per_yr','HeatingWarming_K','Background_C','TotalBasal_C'});

disp(exampleTable);

%% Local function: stable exponential ratio, including x = 0
function y = expRatio(x)
% Returns (exp(x)-1)/x, with its continuous extension y(0)=1.

    y = ones(size(x));
    small = abs(x) < 1e-5;

    xs = x(small);
    y(small) = 1 + xs/2 + xs.^2/6 + xs.^3/24;

    y(~small) = expm1(x(~small))./x(~small);
end
