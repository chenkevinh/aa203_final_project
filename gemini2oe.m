function [a, ecc, incl] = gemini2oe(h_a,h_p,incl,T)

mu = 398600.4418; % km^3/s^2
re = 6378.137; % km

h_a = h_a * 1.852; % km
T = T*60; % s

r_a = h_a + re;
a = ((T/(2*pi))^2*mu)^(1/3);
ecc = r_a/a-1;
incl = deg2rad(incl);

end