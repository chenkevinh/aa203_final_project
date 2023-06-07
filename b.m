function B = b(oe0)

mu = 398600.4418; % km^3/s^2

a = oe0(1);
u = oe0(5)+oe0(6);
n = sqrt(mu/a^3);

B = 1/(a*n)*[0 2 0;
            -2 0 0;
            sin(u) 2*cos(u) 0;
            -cos(u) 2*sin(u) 0;
            0 0 cos(u);
            0 0 sin(u)];

end