clc; clear; close all;

% Simulation

addpath('mean_osc')

% Constants
mu = 398600.4418; % km^3/s^2
re = 6378.137; % km
j2 = 1.08262668e-3;

% GATV initial conditions
a0i = 3603.3; % n. mi.
a0i = a0i * 1.852; % km
ecc0i = 0.0007;
incl0i = 28.87; % deg
incl0i = deg2rad(incl0i);
RAAN0i = 68.2; % deg
RAAN0i = deg2rad(RAAN0i);
argp0i = 94.52; % deg
argp0i = deg2rad(argp0i);
nu0i = 4.34; % deg
nu0i = deg2rad(nu0i);
M0i = true2mean(nu0i,ecc0i);
oe0i_mean = [a0i, ecc0i, incl0i, RAAN0i, argp0i, M0i]';

n = sqrt(mu/a0i^3);

% Gemini Spacecraft initial conditions
h_a = 145.7; % n. mi.
h_p = 144.1; % n. mi.
incl = 28.87;
T = 89.89; % min
[a1i, ecc1i, incl1i] = gemini2oe(h_a,h_p,incl,T);
RAAN1i = RAAN0i;
argp1i = argp0i;
nu1i = nu0i-0.00734;
M1i = true2mean(nu1i,ecc1i);
oe1i_mean = [a1i, ecc1i, incl1i, RAAN1i, argp1i, M1i]';

oe0i_osc = mean2osc(oe0i_mean,1);
oe1i_osc = mean2osc(oe1i_mean,1);

roei = oe2roe(oe0i_mean,oe1i_mean);

[r0i, v0i] = oe2eci(oe0i_osc);
[r1i, v1i] = oe2eci(oe1i_osc);
y0i = [r0i; v0i];
y1i = [r1i; v1i];

T = 2*pi*sqrt(a0i^3/mu); % s
num_orbits = 2;
tstart = 0; tend = num_orbits*T; tint = T/100;
ts = tstart:tint:tend; num_steps = numel(ts);
ts = [ts ts(end) + tint];

y0 = zeros(6,num_steps);
y1 = zeros(6,num_steps);
y0(:,1) = y0i;
y1(:,1) = y1i;
r0_eci = zeros(3,num_steps);
r1_eci = zeros(3,num_steps);
v0_eci = zeros(3,num_steps);
v1_eci = zeros(3,num_steps);

oe0_osc = zeros(6,num_steps);
oe0_mean = oe0_osc;
oe1_osc = oe0_osc;
oe1_mean = oe0_osc;
roe_fode = zeros(6,num_steps);
rho_eci = zeros(3,num_steps);
rhodot_eci = zeros(3,num_steps);
rho_rtn = rho_eci;
rhodot_rtn = rho_rtn;

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
for i = 1:num_steps
    if i ~= 1
        tsteps = ts(i-1:i+1);
        [~, y0_i] = ode113(@(t,y) state_deriv_j2_u(t,y,[0;0;0]),tsteps,y0(:,i-1),options);
        [~, y1_i] = ode113(@(t,y) state_deriv_j2_u(t,y,[0;0;0]),tsteps,y1(:,i-1), options);
        y0(:,i) = y0_i(2,:);
        y1(:,i) = y1_i(2,:);
    end
    r0_eci(:,i) = y0(1:3,i);
    r1_eci(:,i) = y1(1:3,i);
    v0_eci(:,i) = y0(4:6,i);
    v1_eci(:,i) = y1(4:6,i);
    rho_eci(:,i) = r1_eci(:,i) - r0_eci(:,i);
    rhodot_eci(:,i) = v1_eci(:,i) - v0_eci(:,i);
    [rho_rtn(:,i),rhodot_rtn(:,i)] = eci2rtn(rho_eci(:,i),rhodot_eci(:,i),r0_eci(:,i),v0_eci(:,i));
    oe0_osc(:,i) = eci2oe(r0_eci(:,i),v0_eci(:,i));
    oe0_mean(:,i) = osc2mean(oe0_osc(:,i));
    oe1_osc(:,i) = eci2oe(r1_eci(:,i),v1_eci(:,i));
    oe1_mean(:,i) = osc2mean(oe1_osc(:,i));
end
ts = ts(1:end-1);

% Calculate ROE
oe0_temp = oe0_mean;
oe1_temp = oe1_mean;
oe0_temp(6,:) = unwrap(oe0_temp(6,:));
oe1_temp(6,:) = unwrap(oe1_temp(6,:));
for i = 1:num_steps
    roe_fode(:,i) = a0i*1000*oe2roe(oe0_temp(:,i),oe1_temp(:,i));
end

roe_stm = zeros(6,num_steps);

% Calculate ROE from STM
for i = 1:num_steps
    if i == 1
        roe_stm(:,i) = a0i*1000*roei;
    else
        Phi = stm_j2(oe0_mean(:,i-1),tint);
        roe_stm(:,i) = Phi*roe_stm(:,i-1);
    end
end

% Final ROEs
roef = [0; 45; 0; 0; 0; 0];
droe = roef - roei;

% Plot ROEs
labels = {'a\deltaa (m)','a\delta\lambda (m)','a\deltae_x (m)','a\deltae_y (m)','a\deltai_x (m)','a\deltai_y (m)'};
figure()
for i = 1:6
    subplot(3,2,i)
    plot(ts/T,roe_fode(i,:),'LineWidth',2)
    hold on
%     plot(ts/T,roef(i)*ones(1,num_steps),'r--')
    plot(ts/T,roe_stm(i,:),'LineWidth',2)
    ylabel(labels(i))
    xlabel('Orbital Periods')
    grid on
end
set(gcf,'Position',[0 0 1000 350])
% exportgraphics(gcf,'203_figures\roe_stm.eps')

% % Plot STM error
% figure()
% labels = {'a\deltaa Error (m)','a\delta\lambda Error (m)','a\deltae_x Error (m)',...
%           'a\deltae_y Error (m)','a\deltai_x Error (m)','a\deltai_y Error (m)'};
% for i = 1:6
%     subplot(3,2,i)
%     plot(ts/T, roe_stm(i,:)-roe_fode(i,:),'LineWidth',2)
%     grid on
%     ylabel(labels(i))
%     xlabel('Orbital Periods')
% end

%% Plot 3D

[xE, yE, zE] = ellipsoid(0,0,0,re,re,re,20);

figure()
view([175 30])
surface(xE,yE,zE,'FaceColor','#4DBEEE','EdgeColor','black');
axis equal
hold on
plot3(r0_eci(1,:),r0_eci(2,:),r0_eci(3,:),'.r','MarkerSize',20);
plot3(r1_eci(1,:),r1_eci(2,:),r1_eci(3,:),'.b','MarkerSize',20);
xlabel('I (km)')
ylabel('J (km)')
zlabel('K (km)')
grid on
% exportgraphics(gcf,'203_figures\orbit_eci.eps')