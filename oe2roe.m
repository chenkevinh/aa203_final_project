function roe = oe2roe(oe0,oe1)

[a0,ecc0,incl0,RAAN0,argp0,M0] = struct('x',num2cell(oe0)).x;
[a1,ecc1,incl1,RAAN1,argp1,M1] = struct('x',num2cell(oe1)).x;

roe = [(a1-a0)/a0;
        M1 + argp1 - (M0 + argp0) + (RAAN1-RAAN0)*cos(incl0);
        ecc1*cos(argp1)-ecc0*cos(argp0);
        ecc1*sin(argp1)-ecc0*sin(argp0);
        incl1-incl0;
        (RAAN1-RAAN0)*sin(incl0)];

end