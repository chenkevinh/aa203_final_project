function [statedot] = state_deriv_j2(t,state)

% Constants
re = 6378.137; % km
mu = 398600.4418; % km^3/s^2
j2 = 1.08262668e-3;

r_eci = state(1:3);

r = norm(r_eci);

agrav = -mu*r_eci/r^3; % km/s^2

f_j2_eci = mu*r_eci/r^3 .* [3/2*j2*(re/r)^2*(5*r_eci(3)^2/r^2 - 1);
                          3/2*j2*(re/r)^2*(5*r_eci(3)^2/r^2 - 1);
                          3/2*j2*(re/r)^2*(5*r_eci(3)^2/r^2 - 3)];

statedot = zeros(6,1);
statedot(1:3) = state(4:6);
statedot(4:6) = agrav + f_j2_eci;

end