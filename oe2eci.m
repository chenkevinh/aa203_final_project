function [r_eci, v_eci] = oe2eci(oe)

[a,ecc,incl,RAAN,argp,M] = struct('x',num2cell(oe)).x;

nu = mean2true(M,ecc,1e-6);

% if ecc == 0
%     if incl == 0
%         argp = 0;
%         RAAN = 0;
%         nu = truelon;
%     else
%         argp = 0;
%         nu = arglat;
%     end
% end
% 
% if ecc ~= 0 && incl == 0
%     RAAN = 0;
%     argp = lonper;
% end

mu = 398600.4418; % km^3/s^2

p = a*(1-ecc^2);
r = p/(1+ecc*cos(nu));
rPQW = r*[cos(nu) sin(nu) 0]';
R1 = [cos(-argp) sin(-argp) 0;
      -sin(-argp) cos(-argp) 0;
      0 0 1];
R2 = [1 0 0;
      0 cos(-incl) sin(-incl);
      0 -sin(-incl) cos(-incl)];
R3 = [cos(-RAAN) sin(-RAAN) 0;
      -sin(-RAAN) cos(-RAAN) 0;
      0 0 1];
pqw2eci = R3*R2*R1;

r_eci = pqw2eci*rPQW;

vPQW = sqrt(mu/p)*[-sin(nu) ecc+cos(nu) 0]';
v_eci = pqw2eci*vPQW;

end