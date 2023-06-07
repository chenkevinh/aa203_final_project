function [rRTN, vRTN] = eci2rtn(rECI, vECI, r0, v0)

Rhat = r0/norm(r0);
Nhat = cross(r0,v0)/norm(cross(r0,v0));
That = cross(Nhat,Rhat);
R = [Rhat That Nhat]^-1;
rRTN = R*rECI;

h0 = cross(r0,v0);
thetadot = norm(h0)/norm(r0)^2;
vRTN = R*vECI - cross([0 0 thetadot]',rRTN);

end