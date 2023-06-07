clc; clear; close all;

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
nu1i = deg2rad(3.8586);
M1i = true2mean(nu1i,ecc1i);
oe1i_mean = [a1i, ecc1i, incl1i, RAAN1i, argp1i, M1i]';

oe0i_osc = mean2osc(oe0i_mean,1);
oe1i_osc = mean2osc(oe1i_mean,1);

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
y0(:,1) = y0i;
r0_eci = zeros(3,num_steps);
v0_eci = zeros(3,num_steps);
oe0_osc = zeros(6,num_steps);
oe0_mean = oe0_osc;

% FODE simulation of chief orbit with J2
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
for i = 1:num_steps
    if i ~= 1
        tsteps = ts(i-1:i+1);
        [~, y0_i] = ode113(@state_deriv_j2, tsteps, y0(:,i-1), options);
        y0(:,i) = y0_i(2,:);
    end
    r0_eci(:,i) = y0(1:3,i);
    v0_eci(:,i) = y0(4:6,i);
    oe0_osc(:,i) = eci2oe(r0_eci(:,i),v0_eci(:,i));
    oe0_mean(:,i) = osc2mean(oe0_osc(:,i));
end
ts = ts(1:end-1);

% ROEs
roei = oe2roe(oe0i_mean,oe1i_mean);
roef = [0; 45; 0; 0; 0; 0]/a0i/1000;
droe = roef - roei;

% Impulsive maneuvers with Simulation
dvt1 = n*a0i/4*(droe(1) + norm(droe(3:4)));
dvt2 = n*a0i/4*(droe(1) - norm(droe(3:4)));
burnvt1 = dvt1*[0; 1; 0];
burnvt2 = dvt2*[0; 1; 0];
uM1 = wrapTo2Pi(atan(droe(4)/droe(3))) + pi;
uM2 = uM1 + pi;

dvr1 = n*a0i/2*(-droe(2)/2 + norm(droe(3:4)));
dvr2 = n*a0i/2*(-droe(2)/2 - norm(droe(3:4)));
burnvr1 = [dvr1; 0; 0];
burnvr2 = [dvr2; 0; 0];
uM3 = uM2 + pi;
uM4 = uM3 + pi;

us = M0i + argp0i + n*ts;
[~,burn1i] = min(abs(us-uM1));
[~,burn2i] = min(abs(us-uM2));
[~,burn3i] = min(abs(us-uM3));
[~,burn4i] = min(abs(us-uM4));
[~,burnoop] = min(abs(us-(pi/2+3*pi)));

us = zeros(3,num_steps);

% Simulation using STM with maneuvers
roe = zeros(6,num_steps);
for i = 1:num_steps
    if i == 1
        roe(:,i) = roei;
        continue
    end
    
    if i == burn1i
        B = b(oe0_mean(:,i-1));
        u = burnvt1;
    elseif i == burn2i
        B = b(oe0_mean(:,i-1));
        u = burnvt2;
    elseif i == burn3i
        B = b(oe0_mean(:,i-1));
        dvr1 = n*a0i/2*(-(roe(2,i-1)-roef(2))/2 + norm(roe(3:4,i-1)));
        u = [-dvr1; 0; 0];
    elseif i == burn4i
        B = b(oe0_mean(:,i-1));
        dvr2 = n*a0i/2*(-(roe(2,i-1)-roef(2))/2 - norm(roe(3:4,i-1)));
        u = [-dvr2; 0; 0];
    elseif i == burnoop
        B = b(oe0_mean(:,i-1));
        dvoop = -n*a0i*(norm(roe(5:6,i-1)));
        u = [0; 0; dvoop];
    else
        B = zeros(6,3);
        u = zeros(3,1);
    end

    us(:,i) = u;

    Phi = stm_j2(oe0_mean(:,i-1),tint);
    roe(:,i) = Phi*roe(:,i-1) + B*u;
end

% Plot ROEs
labels = {'a\deltaa (m)','a\delta\lambda (m)','a\deltae_x (m)','a\deltae_y (m)','a\deltai_x (m)','a\deltai_y (m)'};
figure()
for i = 1:6
    subplot(3,2,i)
    plot(ts/T,a0i*1000*roe(i,:),'LineWidth',2)
    hold on
    plot(ts/T,a0i*1000*roef(i)*ones(1,num_steps),'r--')
    ylabel(labels(i))
    xlabel('Orbital Periods')
    grid on
end
set(gcf,'Position',[0 0 1000 350])
% exportgraphics(gcf,'203_figures\closed_form.eps')

% Plot deltaV history
figure()
labels = {'{\Delta}v_R (m/s)','{\Delta}v_T (m/s)','{\Delta}v_N (m/s)','{\Delta}v (m/s)'};
for i = 1:3
    subplot(4,2,2*i-1)
    plot(ts(1:end)/T,us(i,:)*1000,'LineWidth',2)
    ylabel(labels(i))
    grid on
end
subplot(4,2,7)
plot(ts(1:end)/T,sum(abs(us),1)*1000,'LineWidth',2)
ylabel(labels(4))
xlabel('Orbital Periods')
grid on

deltav_cum = cumsum(sum(abs(us),1)*1000);

% Plot cumulative deltaV
subplot(4,2,2:2:8)
plot(ts(1:end)/T,deltav_cum,'LineWidth',2)
grid on
ylabel('Cumulative \Deltav (m/s)')
xlabel('Orbital Periods')
set(gcf,'Position',[0 0 1000 350])
% exportgraphics(gcf,'203_figures\closed_form_deltav.eps')

%% iLQR Solution

% Calculate STM
As = zeros(6,6,num_steps);
Bs = zeros(6,3,num_steps);
for i = 1:num_steps
    As(:,:,i) = stm_j2(oe0_mean(:,i),tint);
    Bs(:,:,i) = b(oe0_mean(:,i));
end

n = 6;
m = 3;
Q = 0.1*diag([1, 1, 1, 1, 1, 1]);
R = 10*eye(m);
QN = 1e7*eye(n);
N = num_steps - 1;
[roe_bar, u_bar, Y, y] = ilqr(As, Bs, roei, roef, N, Q, R, QN, 1e-3, 1000);

roe_lqr = zeros(n,N+1);
u = zeros(m,N);
roe_lqr(:,1) = roei;
for k = 1:N
    u(:,k) = u_bar(:,k) + y(:,k) + Y(:,:,k) * (roe_lqr(:,k) - roe_bar(:,k));
    roe_lqr(:,k+1) = As(:,:,k)*roe_lqr(:,k) + Bs(:,:,k)*u(:,k);
end

% Plot ROEs from iLQR solution
labels = {'a\deltaa (m)','a\delta\lambda (m)','a\deltae_x (m)','a\deltae_y (m)','a\deltai_x (m)','a\deltai_y (m)'};
figure()
for i = 1:6
    subplot(3,2,i)
    plot(ts/T,a0i*1000*roe_lqr(i,:),'LineWidth',2)
    hold on
    plot(ts/T,a0i*1000*roef(i)*ones(1,num_steps),'r--')
    ylabel(labels(i))
    xlabel('Orbital Periods')
    grid on
end
set(gcf,'Position',[0 0 1000 400])
% exportgraphics(gcf,'203_figures\roe_ilqr.eps')

% Plot ROEs from iLQR solution
labels = {'a\deltaa (m)','a\delta\lambda (m)','a\deltae_x (m)','a\deltae_y (m)','a\deltai_x (m)','a\deltai_y (m)'};
figure()
for i = 1:6
    subplot(3,2,i)
    plot(ts(end-20:end)/T,a0i*1000*roe_lqr(i,(end-20:end)),'LineWidth',2)
    hold on
    plot(ts(end-20:end)/T,a0i*1000*roef(i)*ones(size(ts(end-20:end))),'r--')
    ylabel(labels(i))
    xlabel('Orbital Periods')
    grid on
end
set(gcf,'Position',[0 0 1000 400])
% exportgraphics(gcf,'203_figures\roe_ilqr_final.eps')

% Plot deltaV history
figure()
labels = {'{\Delta}v_R (m/s)','{\Delta}v_T (m/s)','{\Delta}v_N (m/s)','\Sigma |{\Delta}v| (m/s)'};
for i = 1:3
    subplot(4,2,2*i-1)
    plot(ts(1:end-1)/T,u(i,:)*1000,'LineWidth',2)
    ylabel(labels(i))
    grid on
end
subplot(4,2,7)
plot(ts(1:end-1)/T,sum(abs(u),1)*1000,'LineWidth',2)
ylabel(labels(4))
xlabel('Orbital Periods')
grid on

deltav_cum = cumsum(sum(abs(u),1)*1000);

% Plot cumulative deltaV
subplot(4,2,2:2:8)
plot(ts(1:end-1)/T,deltav_cum,'LineWidth',2)
grid on
ylabel('Cumulative \Deltav (m/s)')
xlabel('Orbital Periods')
set(gcf,'Position',[0 0 1000 350])
% exportgraphics(gcf,'203_figures\deltav_hist.eps')

% deltaV LBs

da_des = droe(1);
dlambda_des = droe(2);
de_des = droe(3:4);
Ms = unwrap(oe0_mean(6,:));
Mf = Ms(end);
Mi = Ms(1);
deltaM = Mf - Mi;

n = sqrt(mu/a0i^3);
deltavLBs = n*a0i*1000*[abs(da_des)/2,abs(dlambda_des)/(3*deltaM),norm(de_des)/2]

% % Add some noise
% roe_lqr = zeros(n,N+1);
% u = zeros(m,N);
% roe_lqr(:,1) = roei;
% for k = 1:N
%     u(:,k) = u_bar(:,k) + y(:,k) + Y(:,:,k) * (roe_lqr(:,k) - roe_bar(:,k));
%     roe_lqr(:,k+1) = As(:,:,k)*roe_lqr(:,k) + Bs(:,:,k)*(u(:,k) + 0.00001*randn(3,1));
% end
% 
% % Plot ROEs
% labels = {'a\deltaa (m)','a\delta\lambda (m)','a\deltae_x (m)','a\deltae_y (m)','a\deltai_x (m)','a\deltai_y (m)'};
% figure()
% for i = 1:6
%     subplot(3,2,i)
%     plot(ts/T,a0i*1000*roe_lqr(i,:),'LineWidth',2)
%     hold on
%     plot(ts/T,roef(i)*ones(1,num_steps),'r--')
%     ylabel(labels(i))
%     xlabel('Orbital Periods')
%     grid on
% end
% 
% % ROE Error
% labels = {'a\deltaa (m)','a\delta\lambda (m)','a\deltae_x (m)','a\deltae_y (m)','a\deltai_x (m)','a\deltai_y (m)'};
% figure()
% for i = 1:6
%     subplot(3,2,i)
%     plot(ts/T,roe_lqr(i,:)-roef(i),'LineWidth',2)
%     ylabel(labels(i))
%     xlabel('Orbital Periods')
%     grid on
% end
% 
% figure()
% labels = {'{\Delta}v_R (m/s)','{\Delta}v_T (m/s)','{\Delta}v_N (m/s)','{\Delta}v (m/s)'};
% for i = 1:3
%     subplot(4,1,i)
%     plot(ts(1:end-1)/T,u(i,:)*1000,'LineWidth',2)
%     ylabel(labels(i))
%     grid on
% end
% subplot(4,1,4)
% plot(ts(1:end-1)/T,vecnorm(u)*1000,'LineWidth',2)
% ylabel(labels(4))
% xlabel('Orbital Periods')
% grid on
% 
% deltav_cum = cumsum(vecnorm(u)*1000);
% 
% figure()
% plot(ts(1:end-1)/T,deltav_cum,'LineWidth',2)
% grid on
% ylabel('Cumulative \Deltav (m/s)')
% xlabel('Orbital Periods')

%% FODE Sim of Maneuvers

y1 = zeros(6,num_steps);
y1(:,1) = y1i;

u_eci = zeros(3,num_steps-1);
for i = 1:num_steps-1
    u_eci(:,i) = rtn2eci(u(:,i),[0;0;0],r0_eci(:,i),v0_eci(:,i));
end

ts = [ts ts(end)+tint];
for i = 2:num_steps
    tsteps = ts(i-1:i+1);
    [~, y1_i] = ode113(@(t,y) state_deriv_j2_u(t,y,u_eci(:,i-1)/tint),tsteps,y1(:,i-1), options);
    y1(:,i) = y1_i(2,:);
end
ts = ts(1:end-1);

r1_eci = y1(1:3,:);
v1_eci = y1(4:6,:);

oe1_osc = zeros(6,num_steps);
oe1_mean = zeros(6,num_steps);
for i = 1:num_steps
    oe1_osc(:,i) = eci2oe(r1_eci(:,i),v1_eci(:,i));
    oe1_mean(:,i) = osc2mean(oe1_osc(:,i));
end

roe_fode = zeros(6,num_steps);
oe0_unwrapped = oe0_mean;
oe0_unwrapped(6,:) = unwrap(oe0_unwrapped(6,:));
oe1_unwrapped = oe1_mean;
oe1_unwrapped(6,:) = unwrap(oe1_unwrapped(6,:));

for i = 1:num_steps
    roe_fode(:,i) = oe2roe(oe0_unwrapped(:,i), oe1_unwrapped(:,i));
end

% Plot ROEs
labels = {'a\deltaa (m)','a\delta\lambda (m)','a\deltae_x (m)','a\deltae_y (m)','a\deltai_x (m)','a\deltai_y (m)'};
figure()
for i = 1:6
    subplot(3,2,i)
    plot(ts/T,a0i*1000*roe_fode(i,:),'LineWidth',2)
    hold on
    plot(ts/T,a0i*1000*roef(i)*ones(1,num_steps),'r--')
    ylabel(labels(i))
    xlabel('Orbital Periods')
    grid on
end

% Calculate Relative Positions and Velocities
rho_eci = r1_eci-r0_eci;
rhodot_eci = v1_eci-v0_eci;
rho_rtn = zeros(3,num_steps);
rhodot_rtn = rho_rtn;
for i = 1:num_steps
    [rho_rtn(:,i),rhodot_rtn(:,i)] = eci2rtn(rho_eci(:,i),rhodot_eci(:,i),r0_eci(:,i),v0_eci(:,i));
end

%% Plot Relative Velocity Time History

figure()
rho = rho_rtn;
labels = {'R (m)','T (m)', 'N (m)'};
for i = 1:3
    subplot(3,1,i)
    plot(ts/T, 1000*rho(i,:),'LineWidth',2)
    grid on
    ylabel(labels(i))
end
% exportgraphics(gcf,'203_figures\velocity_u_maneuver.eps')

%% Plot Relative Positions and Velocities

figure()
rho = rho_rtn;
subplot(2,2,1)
plot3(rho(1,:),rho(2,:),rho(3,:),'LineWidth',2)
xlabel('R (km)')
ylabel('T (km)')
zlabel('N (km)')
axis equal
grid on
subplot(2,2,2)
plot(rho(2,:),rho(1,:),'LineWidth',2)
hold on
xlabel('T (km)')
ylabel('R (km)')
axis equal
grid on
subplot(2,2,3)
plot(rho(3,:),rho(1,:),'LineWidth',2)
hold on
xlabel('N (km)')
ylabel('R (km)')
axis equal
grid on
subplot(2,2,4)
plot(rho(2,:),rho(3,:),'LineWidth',2)
hold on
xlabel('T (km)')
ylabel('N (km)')
axis equal
grid on
% exportgraphics(gcf,'203_figures\position_maneuvers.eps')

figure()
rho = 1000*rhodot_rtn;
subplot(2,2,1)
plot3(rho(1,:),rho(2,:),rho(3,:),'LineWidth',2)
hold on
xlabel('R (m/s)')
ylabel('T (m/s)')
zlabel('N (m/s)')
axis equal
grid on
subplot(2,2,2)
plot(rho(2,:),rho(1,:),'LineWidth',2)
xlabel('T (m/s)')
ylabel('R (m/s)')
axis equal
grid on
subplot(2,2,3)
plot(rho(3,:),rho(1,:),'LineWidth',2)
hold on
xlabel('N (m/s)')
ylabel('R (m/s)')
axis equal
grid on
subplot(2,2,4)
plot(rho(2,:),rho(3,:),'LineWidth',2)
xlabel('T (m/s)')
ylabel('N (m/s)')
axis equal
grid on
% exportgraphics(gcf,'203_figures\velocity_maneuvers.eps')

% %% Plot 3D
% 
% M(num_steps-1) = struct('cdata',[],'colormap',[]);
% [xE, yE, zE] = ellipsoid(0,0,0,re,re,re,20);
% close all
% 
% tic
% parfor i = 1:num_steps-1
%     h = figure();
%     h.Visible = 'off';
%     surface(xE,yE,zE,'FaceColor','#4DBEEE','EdgeColor','black');
%     axis equal
%     hold on
%     plot3(r0_eci(1,i),r0_eci(2,i),r0_eci(3,i),'.r','MarkerSize',20);
%     plot3(r1_eci(1,i),r1_eci(2,i),r1_eci(3,i),'.b','MarkerSize',20);
%     quiver3(r1_eci(1,i),r1_eci(2,i),r1_eci(3,i),u_eci(1,i),u_eci(2,i),u_eci(3,i),10000,'black','LineWidth',2)
%     hold off
%     v = r1_eci(:,i);
%     view(v)
%     camzoom(300)
%     set(gcf,'Position',[0 0 1920 1080])
%     M(i) = getframe;
% end
% toc
% 
% %%
% close all
% figure()
% axis off
% set(gcf,'Position',[0 0 1920 1080])
% movie(M,1,10)