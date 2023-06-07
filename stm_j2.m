function STM = stm_j2(oe0,tau)

mu = 398600.4418; % km^3/s^2
re = 6378.137; % km
j2 = 1.08262668e-3;

[a0,ecc0,incl0,RAAN0i,argp0i,M0] = struct('x',num2cell(oe0)).x;

eta = sqrt(1-ecc0^2);
k = 3/4*j2*re^2*sqrt(mu)/(a0^(7/2)*eta^4);
E = 1+eta;
F = 4+3*eta;
G = 1/eta^2;
P = 3*cos(incl0)^2-1;
Q = 5*cos(incl0)^2-1;
R = cos(incl0);
S = sin(2*incl0);
T = sin(incl0)^2;

argpdot = k*Q;
RAANdot = -2*k*R;

argp0f = argp0i + argpdot*tau;
RAAN0f = RAAN0i + RAANdot*tau;

exi = ecc0*cos(argp0i);
eyi = ecc0*sin(argp0i);
exf = ecc0*cos(argp0f);
eyf = ecc0*sin(argp0f);

n = sqrt(mu/a0^3);

STM = [1 0 0 0 0 0;
           -(3/2*n + 7/2*k*E*P)*tau 1 k*exi*F*G*P*tau k*eyi*F*G*P*tau -k*F*S*tau 0;
           7/2*k*eyf*Q*tau 0 cos(argpdot*tau)-4*k*exi*eyf*G*Q*tau -sin(argpdot*tau)-4*k*eyi*eyf*G*Q*tau 5*k*eyf*S*tau 0;
           -7/2*k*exf*Q*tau 0 sin(argpdot*tau)+4*k*exi*exf*G*Q*tau cos(argpdot*tau)+4*k*eyi*exf*G*Q*tau -5*k*exf*S*tau 0;
           0 0 0 0 1 0;
           7/2*k*S*tau 0 -4*k*exi*G*S*tau -4*k*eyi*G*S*tau 2*k*T*tau 1];

end