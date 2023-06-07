function [rECI, vECI] = rtn2eci(rRTN, vRTN, r0, v0)

Rhat = r0/norm(r0);
Nhat = cross(r0,v0)/norm(cross(r0,v0));
That = cross(Nhat,Rhat);
R = [Rhat That Nhat];
rECI = R*rRTN;

h0 = cross(r0,v0);
thetadot = norm(h0)/norm(r0)^2;
vECI  = R*(vRTN + cross([0 0 thetadot]',rRTN));

end