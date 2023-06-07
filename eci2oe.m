function oe = eci2oe(r_eci, v_eci)

mu = 398600.4418; % km^3/s^2

rmag = norm(r_eci);
vmag = norm(v_eci);

h = cross(r_eci, v_eci); hmag = norm(h);
n = cross([0,0,1]',h);
nmag = norm(n);
evec = ((vmag^2-mu/rmag)*r_eci - dot(r_eci,v_eci)*v_eci)/mu; ecc = norm(evec);

a = (2/rmag - vmag^2/mu)^-1;

incl = acos(h(3)/hmag);

RAAN = acos(n(1)/nmag);
if n(2) < 0
    RAAN = 2*pi-RAAN;
end

% Argument of periapsis (argp/omega)
argp = acos(dot(n,evec)/(nmag*ecc));
if evec(3) < 0
    argp = 2*pi-argp;
end

% True anomaly (nu)
nu = acos(dot(evec,r_eci)/(ecc*rmag));
if dot(r_eci,v_eci) < 0
    nu = 2*pi-nu;
end

% True longitude (l)
if incl == 0 && ecc == 0
    RAAN = 0;
    argp = 0;
    nu = 0;
    truelon = acos(r_eci(1)/rmag);
    if r_eci(2) < 0
        truelon = 2*pi - truelon;
    end
else
    truelon = RAAN + argp + nu;
end

% Argument of latitude (u)
% if ecc == 0 && incl ~= 0
%     argp = 0;
%     nu = 0;
arglat = acos(dot(n,r_eci)/(nmag*rmag));
if r_eci(3) < 0
    arglat = 2*pi - arglat;
end
% else
%     arglat = argp + nu;
% end

% Longitude of Periapsis
if incl == 0 && ecc ~= 0
    RAAN = 0;
    argp = 0;
    lonper = evec(1)/norm(evec);
    if (evec(3) < 0)
        lonper = 2*pi-lonper;
    end
else
    lonper = RAAN + argp;
end

M = true2mean(nu,ecc);

oe = [a; ecc; incl; RAAN; argp; M];

end