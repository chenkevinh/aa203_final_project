
% transformationmatrix_osc2mean_equinoctial Jacobian of nonlinear osc2mean
% function provided by Alfriend.  
%   [ D_J2, osc_c, Psi_Lam ] = transformationmatrix_osc2mean_equinoctial(mean_c) 
%   Provides transformation matrix and mapped osculating equinoctial 
%   orbital elements using the Jacobian of the nonlinear transformation 
%   provided by Alfriend
% 
%   INPUTS:
%   mean_equi: mean equinoctial orbital element vector
% 
%   OUTPUTS:
%   D_J2: Transformation matrix from mean to osculating equioe
%         D_J2 = I + (-J2*Re^2)*(D_lp+D_sp1+D_sp2)
%   osc_equi: Osculating equioe from transformation matrix
%   Psi_Lam: psi-lambda (difference between true long. and mean long.)

function [D_J2, osc_equi, Psi_Lam] = ...
    transformationmatrix_mean2osc_equinoctial(mean_equi)

J2 = 0.0010826358191967;
mu = 3.986004415e14;
Re = 6.378136300e6;

coef = -J2 * Re^2;

a = mean_equi(1);
Psi = mean_equi(2); 
tq1 = mean_equi(3);
tq2 = mean_equi(4);
p1 = mean_equi(5);
p2 = mean_equi(6);

s_Psi = sin(Psi);
c_Psi = cos(Psi);
s_2Psi = sin(2*Psi);
c_2Psi = cos(2*Psi);
s_3Psi = sin(3*Psi);
c_3Psi = cos(3*Psi);
eps1 = sqrt(tq1^2 + tq2^2);
eps2 = tq1*c_Psi + tq2*s_Psi;
eps3 = tq1*sin(Psi) - tq2*c_Psi;
p = a * ( 1 - eps1^2 );
R = p / ( 1 + eps2 );
eta = sqrt( 1 - eps1^2 );
eta_2 = eta^2;
eta_3 = eta^3;
eta_4 = eta_2^2;
eta_6 = eta_3^2;
eta_8 = eta_4^2;
Vr = sqrt(mu/p) * eps3;
Vt = sqrt(mu/p) * ( 1 + eps2 );
sig1 = sqrt( p1^2 + p2^2 );
sig1_2 = sig1^2;
sig1_4 = sig1_2^2;
sig1_6 = sig1^6;
sig1_8 = sig1_4^2;
sig2 = 1 + sig1_2;
sig2_2 = sig2^2;
sig2_3 = sig2^3;
sig3 = 1 - sig1_2;
sig3_2 = sig3^2;

tau1 = p1*c_Psi + p2*s_Psi;
tau2 = p1*s_Psi - p2*c_Psi;
tau12tau22 = tau1^2-tau2^2;
xi1 = p1*tq1 + p2*tq2;
xi2 = p1*tq2 - p2*tq1;
xj1 = p1*tq1 - p2*tq2;
xj2 = p1*tq2 + p2*tq1;
xk1 = tq1*xi1 - tq2*xi2;
xk2 = tq1*xi2 + tq2*xi1;
sig32sig12 = sig3_2 - sig1_2;
sig322sig12 = sig3_2-2*sig1_2;
a2eta4sig22 = a^2*eta_4*sig2_2;
a2eta4sig23 = a^2*eta_4*sig2_3;
a2eta4sig2 = a^2*eta_4*sig2;
a2eta6sig22 = a^2*eta_6*sig2_2;
a2eta6sig2 = a^2*eta_6*sig2;

TTheta = 1 + 5*sig3_2/(2*sig32sig12);
PPhi = 28 - 150*sig1_2 + 290*sig1_4 - 215*sig1_6 + 60*sig1_8 - 7*sig1^10;

% Equation of Center
Lambda = truelongitude2meanlongitude(a, Psi, tq1, tq2);
Psi_Lam = Psi - Lambda;
Lam_tq1 = (tq1*Vr)/(eta*Vt) + tq2/(eta*(1+eta)) - eta*R*(a+R)*(tq2+s_Psi)/(p^2);
Lam_tq2 = (tq2*Vr)/(eta*Vt) - tq1/(eta*(1+eta)) + eta*R*(a+R)*(tq1+c_Psi)/(p^2);

% Identity Matrix
DI = eye(6);

% Long period part,  D_lp
Lam_lp = ( xi1*xi2/(8*a^2*eta_4*(1+eta)*sig2_2*sig32sig12^2) )*( 4*eta_2*sig32sig12^2*TTheta ...
   +(1+eta)*PPhi );

a_lp = 0;
Psi_lp = Lam_lp - ( TTheta/(4*a^2*eta_4*(1+eta)*sig2_2) )*( (3*(1+eta)+2*eta_2)*xi1*xi2 ...
   +(1+eta)*(2*(tau1*xi2+tau2*xi1)+eps1^2*tau1*tau2) );
tq1_lp = -( 1/(8*a2eta4sig22*sig32sig12^2) )*( 2*eta_2*(p1*xi1+p2*xi2)*sig32sig12^2*TTheta ...
   +tq2*xi1*xi2*PPhi );
tq2_lp = ( 1/(8*a2eta4sig22*sig32sig12^2) )*( 2*eta_2*(p1*xi2-p2*xi1)*sig32sig12^2*TTheta ...
   +tq1*xi1*xi2*PPhi );
p1_lp = -( sig3/(16*a2eta4sig2*sig32sig12^2) )*( 5*p2*sig2_2*xi1*xi2 ...
   -xk1*sig32sig12^2*TTheta );
p2_lp = ( sig3/(16*a2eta4sig2*sig32sig12^2) )*( 5*p1*sig2_2*xi1*xi2 ...
   +xk2*sig32sig12^2*TTheta );

D_lp_11 = -(1/a) * a_lp;
D_lp_12 = 0;
D_lp_13 = 0;
D_lp_14 = 0;
D_lp_15 = 0;
D_lp_16 = 0;

D_lp_21 = -(2/a) * Psi_lp;
D_lp_22 = -( TTheta/(4*a2eta4sig22) )*( 2*(tau1*xi1-tau2*xi2)+eps1^2*tau12tau22 );
D_lp_23 = -( 1/(8*a2eta6sig22*sig32sig12^2) )*( (4*tq1*xi1*xi2 ...
   +eta_2*(p1*xi2-p2*xi1))*(6*sig32sig12^2*TTheta-PPhi) ...
   +4*(eta_2*(p1*tau2-p2*tau1)+4*tq1*(tau1*xi2+tau2*xi1)+tq1*(2-eta_2)*tau1*tau2)*sig32sig12^2*TTheta );
D_lp_24 = -( 1/(8*a2eta6sig22*sig32sig12^2) )*( (4*tq2*xi1*xi2 ...
   +eta_2*(p1*xi1+p2*xi2))*(6*sig32sig12^2*TTheta-PPhi) ...
   +4*(eta_2*(p1*tau1+p2*tau2)+4*tq2*(tau1*xi2+tau2*xi1)+tq2*(2-eta_2)*tau1*tau2)*sig32sig12^2*TTheta );
D_lp_25 = -( 1/(8*a2eta4sig23*sig32sig12^3) )*( 2*p1*sig32sig12*(5*sig2_2*sig3 ...
   -4*sig32sig12^2*TTheta)*(3*xi1*xi2+2*(tau1*xi2+tau2*xi1)+eps1^2*tau1*tau2) ...
   +2*sig2*sig32sig12^3*TTheta*(3*xk2+4*(xi1*s_Psi+xi2*c_Psi) ...
   +eps1^2*(tau1*s_Psi+tau2*c_Psi)) -sig2*sig32sig12*(xk2*PPhi ...
   -10*p1*xi1*xi2*(30-116*sig1_2+129*sig1_4-48*sig1_6+7*sig1_8))-4*p1*xi1*xi2*(2+sig1_2*(1+3*sig3))*PPhi );
D_lp_26 = -( 1/(8*a2eta4sig23*sig32sig12^3) )*( 2*p2*sig32sig12*(5*sig2_2*sig3 ...
   -4*sig32sig12^2*TTheta)*(3*xi1*xi2+2*(tau1*xi2+tau2*xi1)+eps1^2*tau1*tau2) ...
   -2*sig2*sig32sig12^3*TTheta*(3*xk1+4*(xi1*c_Psi-xi2*s_Psi) ...
   +eps1^2*(tau1*c_Psi-tau2*s_Psi)) +sig2*sig32sig12*(xk1*PPhi ...
   +10*p2*xi1*xi2*(30-116*sig1_2+129*sig1_4-48*sig1_6+7*sig1_8))-4*p2*xi1*xi2*(2+sig1_2*(1+3*sig3))*PPhi );

D_lp_31 = -(2/a) * tq1_lp;
D_lp_32 = 0;
D_lp_33 = -( 1/(8*a2eta6sig22*sig32sig12^2) )*( 2*eta_2*(2*tq1*(p1*xi1+p2*xi2) ...
   +eta_2*(p1^2-p2^2))*sig32sig12^2*TTheta+tq2*(4*tq1*xi1*xi2+eta_2*(p1*xi2-p2*xi1))*PPhi );
D_lp_34 = -( 1/(8*a2eta6sig22*sig32sig12^2) )*( 4*eta_2*(tq2*(p1*xi1+p2*xi2) ...
   +eta_2*p1*p2)*sig32sig12^2*TTheta+((eta_2+4*tq2^2)*xi1*xi2+eta_2*tq2*(p1*xi1+p2*xi2))*PPhi );
D_lp_35 = -( 5*p1/(4*a2eta4sig22*sig32sig12^2) )*( eta_2*sig2*sig3*(p1*xi1+p2*xi2)  ...
   -tq2*xi1*xi2*(30-116*sig1_2+129*sig1_4-48*sig1_6+7*sig1_8) ) ...
   +( TTheta/(2*a^2*eta_2*sig2_3) )*( 2*p1*(p1*xi1+p2*xi2)-sig2*xi1 ) ...
   -( tq2*PPhi/(8*a2eta4sig23*sig32sig12^3) )*( sig2*sig32sig12*xk2 ...
   +4*p1*xi1*xi2*(2+sig1_2*(1+3*sig3)) );
D_lp_36 = -( 5*p2/(4*a2eta4sig22*sig32sig12^2) )*( eta_2*sig2*sig3*(p1*xi1+p2*xi2)  ...
   -tq2*xi1*xi2*(30-116*sig1_2+129*sig1_4-48*sig1_6+7*sig1_8) ) ...
   +( TTheta/(2*a^2*eta_2*sig2_3) )*( 2*p2*(p1*xi1+p2*xi2)-sig2*xi2 ) ...
   +( tq2*PPhi/(8*a2eta4sig23*sig32sig12^3) )*( sig2*sig32sig12*xk1 ...
   -4*p2*xi1*xi2*(2+sig1_2*(1+3*sig3)) );

D_lp_41 = -(2/a) * tq2_lp;
D_lp_42 = 0;
D_lp_43 = ( 1/(8*a2eta6sig22*sig32sig12^2) )*( 4*eta_2*(tq1*(p1*xi2-p2*xi1) ...
   -eta_2*p1*p2)*sig32sig12^2*TTheta+((eta_2+4*tq1^2)*xi1*xi2+eta_2*tq1*(p1*xi2-p2*xi1))*PPhi );
D_lp_44 = ( 1/(8*a2eta6sig22*sig32sig12^2) )*( 2*eta_2*(2*tq2*(p1*xi2-p2*xi1) ...
   +eta_2*(p1^2-p2^2))*sig32sig12^2*TTheta+tq1*(4*tq2*xi1*xi2+eta_2*(p1*xi1+p2*xi2))*PPhi );
D_lp_45 = ( 5*p1/(4*a2eta4sig22*sig32sig12^2) )*( eta_2*sig2*sig3*(p1*xi2-p2*xi1)  ...
   -tq1*xi1*xi2*(30-116*sig1_2+129*sig1_4-48*sig1_6+7*sig1_8) ) ...
   -( TTheta/(2*a^2*eta_2*sig2_3) )*( 2*p1*(p1*xi2-p2*xi1)-sig2*xi2 ) ...
   +( tq1*PPhi/(8*a2eta4sig23*sig32sig12^3) )*( sig2*sig32sig12*xk2 ...
   +4*p1*xi1*xi2*(2+sig1_2*(1+3*sig3)) );
D_lp_46 = ( 5*p2/(4*a2eta4sig22*sig32sig12^2) )*( eta_2*sig2*sig3*(p1*xi2-p2*xi1)  ...
   -tq1*xi1*xi2*(30-116*sig1_2+129*sig1_4-48*sig1_6+7*sig1_8) ) ...
   -( TTheta/(2*a^2*eta_2*sig2_3) )*( 2*p2*(p1*xi2-p2*xi1)+sig2*xi1 ) ...
   -( tq1*PPhi/(8*a2eta4sig23*sig32sig12^3) )*( sig2*sig32sig12*xk1 ...
   -4*p2*xi1*xi2*(2+sig1_2*(1+3*sig3)) );

D_lp_51 = -(2/a) * p1_lp;
D_lp_52 = 0;
D_lp_53 = -( sig3/(16*a2eta6sig2*sig32sig12^2) )*( 5*p2*sig2_2*(4*tq1*xi1*xi2+eta_2*(p1*xi2-p2*xi1)) ...
   -2*(2*tq1*xk1+eta_2*xi1)*sig32sig12^2*TTheta );
D_lp_54 = -( sig3/(16*a2eta6sig2*sig32sig12^2) )*( 5*p2*sig2_2*(4*tq2*xi1*xi2+eta_2*(p1*xi1+p2*xi2)) ...
   -2*(2*tq2*xk1-eta_2*xi2)*sig32sig12^2*TTheta );
D_lp_55 = -( 1/(16*a2eta4sig22*sig32sig12^3) )*( 20*p1*p2*sig2_2*xi1*xi2*(1+sig3_2*(1+sig2)) ...
   -5*sig2_2*sig3*sig32sig12*(p1*sig3*xk1-p2*sig2*xk2) ...
   +(4*p1*xk1-sig2*sig3*(tq1^2-tq2^2))*sig32sig12^3*TTheta );
D_lp_56 = -( 1/(16*a2eta4sig22*sig32sig12^3) )*( 20*p2^2*sig2_2*xi1*xi2*(1+sig3_2*(1+sig2)) ...
   -5*sig2_2*sig3*sig32sig12*(2*p2*xk1-sig2*xi1*xi2) ...
   +2*(2*p2*xk1-tq1*tq2*sig2*sig3)*sig32sig12^3*TTheta );

D_lp_61 = -(2/a) * p2_lp;
D_lp_62 = 0;
D_lp_63 = ( sig3/(16*a2eta6sig2*sig32sig12^2) )*( 5*p1*sig2_2*(4*tq1*xi1*xi2+eta_2*(p1*xi2-p2*xi1)) ...
   +2*(2*tq1*xk2+eta_2*xi2)*sig32sig12^2*TTheta );
D_lp_64 = ( sig3/(16*a2eta6sig2*sig32sig12^2) )*( 5*p1*sig2_2*(4*tq2*xi1*xi2+eta_2*(p1*xi1+p2*xi2)) ...
   +2*(2*tq2*xk2+eta_2*xi1)*sig32sig12^2*TTheta );
D_lp_65 = ( 1/(16*a2eta4sig22*sig32sig12^3) )*( 20*p1^2*sig2_2*xi1*xi2*(1+sig3_2*(1+sig2)) ...
   +5*sig2_2*sig3*sig32sig12*(2*p1*xk2+sig2*xi1*xi2) ...
   -2*(2*p1*xk2-tq1*tq2*sig2*sig3)*sig32sig12^3*TTheta );
D_lp_66 = ( 1/(16*a2eta4sig22*sig32sig12^3) )*( 20*p1*p2*sig2_2*xi1*xi2*(1+sig3_2*(1+sig2)) ...
   -5*sig2_2*sig3*sig32sig12*(p1*sig2*xk1-p2*sig3*xk2) ...
   -(4*p2*xk2+sig2*sig3*(tq1^2-tq2^2))*sig32sig12^3*TTheta );

D_lp = [ D_lp_11  D_lp_12  D_lp_13  D_lp_14  D_lp_15  D_lp_16;
         D_lp_21  D_lp_22  D_lp_23  D_lp_24  D_lp_25  D_lp_26;
         D_lp_31  D_lp_32  D_lp_33  D_lp_34  D_lp_35  D_lp_36;
         D_lp_41  D_lp_42  D_lp_43  D_lp_44  D_lp_45  D_lp_46;
         D_lp_51  D_lp_52  D_lp_53  D_lp_54  D_lp_55  D_lp_56;
         D_lp_61  D_lp_62  D_lp_63  D_lp_64  D_lp_65  D_lp_66 ];
       
       
% First short period part,  D_sp1
Lam_sp1 = -( eps3*sig322sig12/(2*a^2*eta_4*(1+eta)*sig2_2) )*( (1+eps2)*(2+eps2)+eta_2 ) ...
   - ( 3*(3*sig3_2-2)/(2*a2eta4sig22) )*( Psi_Lam+eps3 );

a_sp1 = -( sig322sig12/(a*eta_6*sig2_2) )*( (1+eps2)^3-eta_3 );
Psi_sp1 = Lam_sp1 + ( eps3*sig322sig12/(2*a^2*eta_4*(1+eta)*sig2_2) )*( (1+eps2)^2+eta*(1+eta) );
tq1_sp1 = -( tq1*sig322sig12/(2*a^2*eta_2*(1+eta)*sig2_2) ) ...
   -( sig322sig12/(2*a2eta4sig22) )*( tq1*(1+eps2)+(eta_2+(1+eps2)*(2+eps2))*c_Psi ) ...
   +( 3*tq2*(3*sig3_2-2)/(2*a2eta4sig22) )*( Psi_Lam+eps3 );
tq2_sp1 = -( tq2*sig322sig12/(2*a^2*eta_2*(1+eta)*sig2_2) ) ...
   -( sig322sig12/(2*a2eta4sig22) )*( tq2*(1+eps2)+(eta_2+(1+eps2)*(2+eps2))*s_Psi ) ...
   -( 3*tq1*(3*sig3_2-2)/(2*a2eta4sig22) )*( Psi_Lam+eps3 );
p1_sp1 = -( 3*p2*sig3/(2*a2eta4sig2) )*( Psi_Lam+eps3 );
p2_sp1 =  ( 3*p1*sig3/(2*a2eta4sig2) )*( Psi_Lam+eps3 );

D_sp1_11 = -(1/a) * a_sp1;
D_sp1_12 = ( 3*sig322sig12/(a*eta_6*sig2_2) )*( eps3*(1+eps2)^2 );
D_sp1_13 = -( 3*sig322sig12/(a*eta_8*sig2_2) )*( (1+eps2)^2*(2*tq1*(1+eps2)+eta_2*c_Psi)-tq1*eta_3 );
D_sp1_14 = -( 3*sig322sig12/(a*eta_8*sig2_2) )*( (1+eps2)^2*(2*tq2*(1+eps2)+eta_2*s_Psi)-tq2*eta_3 );
D_sp1_15 = ( 12*p1*sig3/(a*eta_6*sig2_3) )*( (1+eps2)^3-eta_3 );
D_sp1_16 = ( 12*p2*sig3/(a*eta_6*sig2_3) )*( (1+eps2)^3-eta_3 );

D_sp1_21 = -(2/a) * Psi_sp1;
D_sp1_22 = ( sig322sig12/(2*a^2*eta_4*(1+eta)*sig2_2) )*( eps2*(eta-(1+eps2))+eps3^2 ) ...
   - ( 3*(3*sig3_2-2)/(2*a2eta4sig22*(1+eps2)^2) )*( (1+eps2)^3-eta_3 );
D_sp1_23 = -( 3*(3*sig3_2-2)/(2*a2eta6sig22) )*( 4*tq1*(Psi_Lam+eps3)+eta_2*s_Psi ) ...
   +( 3*(3*sig3_2-2)/(2*a2eta4sig22) )*( Lam_tq1 ) ...
   +( sig322sig12/(2*a^2*eta_6*(1+eta)^2*sig2_2) )*( (eta-(1+eps2))*(tq1*eps3*(4+5*eta) ...
   +eta_2*(1+eta)*s_Psi)-eta*(1+eta)*eps3*(tq1+eta*c_Psi) );
D_sp1_24 = -( 3*(3*sig3_2-2)/(2*a2eta6sig22) )*( 4*tq2*(Psi_Lam+eps3)-eta_2*c_Psi ) ...
   +( 3*(3*sig3_2-2)/(2*a2eta4sig22) )*( Lam_tq2 ) ...
   +( sig322sig12/(2*a^2*eta_6*(1+eta)^2*sig2_2) )*( (eta-(1+eps2))*(tq2*eps3*(4+5*eta) ...
   -eta_2*(1+eta)*c_Psi)-eta*(1+eta)*eps3*(tq2+eta*s_Psi) );
D_sp1_25 = -( 6*p1*sig3*eps3/(a^2*eta_4*(1+eta)*sig2_3) )*( eta-(1+eps2) ) ...
   -( 12*p1*(1-3*sig3)/(a2eta4sig23) )*( Psi_Lam+eps3 );
D_sp1_26 = -( 6*p2*sig3*eps3/(a^2*eta_4*(1+eta)*sig2_3) )*( eta-(1+eps2) ) ...
   -( 12*p2*(1-3*sig3)/(a2eta4sig23) )*( Psi_Lam+eps3 );

D_sp1_31 = -(2/a) * tq1_sp1;
D_sp1_32 = ( sig322sig12/(2*a2eta4sig22) )*( eps3*(tq1+(3+2*eps2)*c_Psi) ...
   +(eta_2+(1+eps2)*(2+eps2))*s_Psi ) ...
   +( 3*tq2*(3*sig3_2-2)/(2*a2eta4sig22*(1+eps2)^2) )*( (1+eps2)^3-eta_3 );
D_sp1_33 = ( 3*tq2*(3*sig3_2-2)/(2*a2eta6sig22) )*( 4*tq1*(Psi_Lam+eps3)+eta_2*s_Psi ) ...
   -( 3*tq2*(3*sig3_2-2)/(2*a2eta4sig22) )*( Lam_tq1 ) ...
   -( sig322sig12/(2*a^2*eta_4*(1+eta)^2*sig2_2) )*( tq1^2*(2+3*eta)+eta_2*(1+eta) ) ...
   -( sig322sig12/(2*a2eta6sig22) )*( (1+eps2)*(4*tq1^2+4*tq1*c_Psi*(2+eps2) ...
   +eta_2*(1+2*c_Psi^2))+eta_2*(3*tq1+c_Psi)*c_Psi );
D_sp1_34 = ( 3*(3*sig3_2-2)/(2*a2eta6sig22) )*( (eta_2+4*tq2^2)*(Psi_Lam+eps3)-eta_2*tq2*c_Psi ) ...
   -( 3*tq2*(3*sig3_2-2)/(2*a2eta4sig22) )*( Lam_tq2 ) ...
   -( tq1*tq2*(2+3*eta)*sig322sig12/(2*a^2*eta_4*(1+eta)^2*sig2_2) ) ...
   -( sig322sig12/(2*a2eta6sig22) )*( 2*(1+eps2)*(2*tq1*tq2+2*tq2*c_Psi*(2+eps2) ...
   +eta_2*s_Psi*c_Psi)+eta_2*((tq1+c_Psi)*s_Psi+2*tq2*c_Psi) );
D_sp1_35 = ( 6*tq1*p1*sig3/(a^2*eta_2*(1+eta)*sig2_3) ) +( 6*p1*sig3/(a2eta4sig23) )*( tq1*(1+eps2) ...
   +(eta_2+(1+eps2)*(2+eps2))*c_Psi ) -( 12*tq2*p1*(2-3*sig1_2)/(a2eta4sig23) )*( Psi_Lam+eps3 );
D_sp1_36 = ( 6*tq1*p2*sig3/(a^2*eta_2*(1+eta)*sig2_3) ) +( 6*p2*sig3/(a2eta4sig23) )*( tq1*(1+eps2) ...
   +(eta_2+(1+eps2)*(2+eps2))*c_Psi ) -( 12*tq2*p2*(2-3*sig1_2)/(a2eta4sig23) )*( Psi_Lam+eps3 );

D_sp1_41 = -(2/a) * tq2_sp1;
D_sp1_42 = ( sig322sig12/(2*a2eta4sig22) )*( eps3*(tq2+(3+2*eps2)*s_Psi) ...
   -(eta_2+(1+eps2)*(2+eps2))*c_Psi ) ...
   -( 3*tq1*(3*sig3_2-2)/(2*a2eta4sig22*(1+eps2)^2) )*( (1+eps2)^3-eta_3 );
D_sp1_43 = -( 3*(3*sig3_2-2)/(2*a2eta6sig22) )*( (eta_2+4*tq1^2)*(Psi_Lam+eps3)+eta_2*tq1*s_Psi ) ...
   +( 3*tq1*(3*sig3_2-2)/(2*a2eta4sig22) )*( Lam_tq1 ) ...
   -( tq1*tq2*(2+3*eta)*sig322sig12/(2*a^2*eta_4*(1+eta)^2*sig2_2) ) ...
   -( sig322sig12/(2*a2eta6sig22) )*( 2*(1+eps2)*(2*tq1*tq2+2*tq1*s_Psi*(2+eps2) ...
   +eta_2*s_Psi*c_Psi)+eta_2*(2*tq1*s_Psi+(tq2+s_Psi)*c_Psi) );
D_sp1_44 = -( 3*tq1*(3*sig3_2-2)/(2*a2eta6sig22) )*( 4*tq2*(Psi_Lam+eps3)-eta_2*c_Psi ) ...
   +( 3*tq1*(3*sig3_2-2)/(2*a2eta4sig22) )*( Lam_tq2 ) ...
   -( sig322sig12/(2*a^2*eta_4*(1+eta)^2*sig2_2) )*( tq2^2*(2+3*eta)+eta_2*(1+eta) ) ...
   -( sig322sig12/(2*a2eta6sig22) )*( (1+eps2)*(4*tq2^2+4*tq2*s_Psi*(2+eps2) ...
   +eta_2*(1+2*s_Psi^2))+eta_2*(3*tq2+s_Psi)*s_Psi );
D_sp1_45 = ( 6*tq2*p1*sig3/(a^2*eta_2*(1+eta)*sig2_3) ) +( 6*p1*sig3/(a2eta4sig23) )*( tq2*(1+eps2) ...
   +(eta_2+(1+eps2)*(2+eps2))*s_Psi ) +( 12*tq1*p1*(2-3*sig1_2)/(a2eta4sig23) )*( Psi_Lam+eps3 );
D_sp1_46 = ( 6*tq2*p2*sig3/(a^2*eta_2*(1+eta)*sig2_3) ) +( 6*p2*sig3/(a2eta4sig23) )*( tq2*(1+eps2) ...
   +(eta_2+(1+eps2)*(2+eps2))*s_Psi ) +( 12*tq1*p2*(2-3*sig1_2)/(a2eta4sig23) )*( Psi_Lam+eps3 );

D_sp1_51 = -(2/a) * p1_sp1;
D_sp1_52 = -( 3*p2*sig3/(2*a2eta4sig2*(1+eps2)^2) )*( (1+eps2)^3-eta_3 );
D_sp1_53 = -( 3*p2*sig3/(2*a2eta6sig2) )*( 4*tq1*(Psi_Lam+eps3)+eta_2*s_Psi ) ...
   +( 3*p2*sig3/(2*a2eta4sig2) )*( Lam_tq1 );
D_sp1_54 = -( 3*p2*sig3/(2*a2eta6sig2) )*( 4*tq2*(Psi_Lam+eps3)-eta_2*c_Psi ) ...
   +( 3*p2*sig3/(2*a2eta4sig2) )*( Lam_tq2 );
D_sp1_55 = ( 6*p1*p2/(a2eta4sig22) )*( Psi_Lam+eps3 );
D_sp1_56 = ( 3/(2*a2eta4sig22) )*(4*p2^2-sig2*sig3)*( Psi_Lam+eps3 );

D_sp1_61 = -(2/a) * p2_sp1;
D_sp1_62 = ( 3*p1*sig3/(2*a2eta4sig2*(1+eps2)^2) )*( (1+eps2)^3-eta_3 );
D_sp1_63 = ( 3*p1*sig3/(2*a2eta6sig2) )*( 4*tq1*(Psi_Lam+eps3)+eta_2*s_Psi ) ...
   -( 3*p1*sig3/(2*a2eta4sig2) )*( Lam_tq1 );
D_sp1_64 = ( 3*p1*sig3/(2*a2eta6sig2) )*( 4*tq2*(Psi_Lam+eps3)-eta_2*c_Psi ) ...
   -( 3*p1*sig3/(2*a2eta4sig2) )*( Lam_tq2 );
D_sp1_65 = -( 3/(2*a2eta4sig22) )*(4*p1^2-sig2*sig3)*( Psi_Lam+eps3 );
D_sp1_66 = -( 6*p1*p2/(a2eta4sig22) )*( Psi_Lam+eps3 );

D_sp1 = [ D_sp1_11  D_sp1_12  D_sp1_13  D_sp1_14  D_sp1_15  D_sp1_16;
          D_sp1_21  D_sp1_22  D_sp1_23  D_sp1_24  D_sp1_25  D_sp1_26;
          D_sp1_31  D_sp1_32  D_sp1_33  D_sp1_34  D_sp1_35  D_sp1_36;
          D_sp1_41  D_sp1_42  D_sp1_43  D_sp1_44  D_sp1_45  D_sp1_46;
          D_sp1_51  D_sp1_52  D_sp1_53  D_sp1_54  D_sp1_55  D_sp1_56;
          D_sp1_61  D_sp1_62  D_sp1_63  D_sp1_64  D_sp1_65  D_sp1_66 ];
       
       
% Second short period part,  D_sp2
Lam_sp2 = -( 1/(2*a^2*eta_4*(1+eta)*sig2_2) )*( 6*(eps3*(1+eps2)*(2+eps2)*tau12tau22 ...
   +(1+eta)*(4-sig1_2)*tau1*tau2)+(eta_2+(1+eta)*(4-sig1_2))*(3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) );

a_sp2 = -( 6*(1+eps2)^3/(a*eta_6*sig2_2) )*( tau1^2-tau2^2 );
Psi_sp2 = Lam_sp2 + ( 1/(2*a^2*eta_4*(1+eta)*sig2_2) )*( 6*eps3*(1+eps2)*(2+eps2)*tau12tau22 ...
   +8*(1+eta)*(2+eps2)*tau1*tau2+(1+eta*(1+eta))*(3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) );
tq1_sp2 = ( 1/(4*a2eta4sig22) )*( 2*tq2*(4-sig1_2)*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   -3*tau12tau22*(10*tq1+3*(tq1*eps2+tq2*eps3)+6*(tq1*c_2Psi+tq2*s_2Psi) ...
   +(tq1^2-tq2^2)*c_3Psi+2*tq1*tq2*s_3Psi) ...
   -6*(3-2*eta_2)*(p1*tau1+p2*tau2)-2*(9-2*eta_2)*((p1^2-p2^2)*c_3Psi+2*p1*p2*s_3Psi) );
tq2_sp2 = -( 1/(4*a2eta4sig22) )*( 2*tq1*(4-sig1_2)*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   +3*tau12tau22*(10*tq2-3*(tq1*eps3-tq2*eps2)+6*(tq1*s_2Psi-tq2*c_2Psi) ...
   +(tq1^2-tq2^2)*s_3Psi-2*tq1*tq2*c_3Psi) ...
   -6*(3-2*eta_2)*(p1*tau2-p2*tau1)+2*(9-2*eta_2)*((p1^2-p2^2)*s_3Psi-2*p1*p2*c_3Psi) );
p1_sp2 = -( sig3/(4*a2eta4sig2) )*( 3*(xi1*c_Psi-xi2*s_Psi)+3*(p1*c_2Psi+p2*s_2Psi) ...
   +xj1*c_3Psi+xj2*s_3Psi );
p2_sp2 =  -( sig3/(4*a2eta4sig2) )*( 3*(xi1*s_Psi+xi2*c_Psi)+3*(p1*s_2Psi-p2*c_2Psi) ...
   +xj1*s_3Psi-xj2*c_3Psi );

D_sp2_11 = -(1/a)*a_sp2;
D_sp2_12 = ( 6*(1+eps2)^2/(a*eta_6*sig2_2) )*( 3*eps3*tau12tau22+4*tau1*tau2*(1+eps2) );
D_sp2_13 = -( 18*(1+eps2)^2/(a*eta_8*sig2_2) )*tau12tau22*( 2*tq1*(1+eps2)+eta_2*c_Psi );
D_sp2_14 = -( 18*(1+eps2)^2/(a*eta_8*sig2_2) )*tau12tau22*( 2*tq2*(1+eps2)+eta_2*s_Psi );
D_sp2_15 = ( 12*(1+eps2)^3/(a*eta_6*sig2_3) )*( 2*p1*tau12tau22-sig2*(tau1*c_Psi-tau2*s_Psi) );
D_sp2_16 = ( 12*(1+eps2)^3/(a*eta_6*sig2_3) )*( 2*p2*tau12tau22-sig2*(tau1*s_Psi+tau2*c_Psi) );

D_sp2_21 = -(2/a) * Psi_sp2;
D_sp2_22 = ( 1/(2*a2eta4sig22) )*( 2*tau12tau22*(4*(1+eps2)-(5+3*sig3))-8*eps3*tau1*tau2 ...
   -3*(2+sig3)*((tau1*xi1-tau2*xi2)+(p1^2-p2^2)*(tq1*c_3Psi+tq2*s_3Psi) ...
   +2*p1*p2*(tq1*s_3Psi-tq2*c_3Psi)) );
D_sp2_23 = ( 1/(2*a2eta6sig22) )*( 8*(4*tq1*(1+eps2)-tq1*(5+3*sig3)+eta_2*c_Psi)*tau1*tau2 ...
   -4*tq1*(2+sig3)*(3*(tau1*xi2+tau2*xi1)+(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi) ...
   -2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   -eta_2*(2+sig3)*(3*(p1*tau2-p2*tau1)+(p1^2-p2^2)*s_3Psi-2*p1*p2*c_3Psi) );
D_sp2_24 = ( 1/(2*a2eta6sig22) )*( 8*(4*tq2*(1+eps2)-tq2*(5+3*sig3)+eta_2*s_Psi)*tau1*tau2 ...
   -4*tq2*(2+sig3)*(3*(tau1*xi2+tau2*xi1)+(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi) ...
   -2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   -eta_2*(2+sig3)*(3*(p1*tau1+p2*tau2)-(p1^2-p2^2)*c_3Psi-2*p1*p2*s_3Psi) );
D_sp2_25 = -( p1/(a2eta4sig23) )*( 2*(8*(1+eps2)-(16+3*sig3))*tau1*tau2-(7-sig1_2)*(3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ) ...
   +( 1/(2*a2eta4sig22) )*( 2*(4*(1+eps2)-(5+3*sig3))*(tau1*s_Psi+tau2*c_Psi) ...
   -(2+sig3)*(3*(tq1*tau2+tq2*tau1)+3*(xi1*s_Psi+xi2*c_Psi)+2*xj1*s_3Psi ...
   -2*xj2*c_3Psi) );
D_sp2_26 = -( p2/(a2eta4sig23) )*( 2*(8*(1+eps2)-(16+3*sig3))*tau1*tau2-(7-sig1_2)*(3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ) ...
   -( 1/(2*a2eta4sig22) )*( 2*(4*(1+eps2)-(5+3*sig3))*(tau1*c_Psi-tau2*s_Psi) ...
   -(2+sig3)*(3*(tq1*tau1-tq2*tau2)+3*(xi1*c_Psi-xi2*s_Psi)+2*xj1*c_3Psi ...
   +2*xj2*s_3Psi) );

D_sp2_31 = -(2/a) * tq1_sp2;
D_sp2_32 = ( 3/(4*a2eta4sig22) )*( 2*tq2*(4-sig1_2)*(2*tau12tau22+(tau1*xi1-tau2*xi2) ...
   +(p1^2-p2^2)*(tq1*c_3Psi+tq2*s_3Psi)+2*p1*p2*(tq1*s_3Psi-tq2*c_3Psi)) ...
   +4*tau1*tau2*(10*tq1+3*(tq1*eps2+tq2*eps3)+6*(tq1*c_2Psi+tq2*s_2Psi)+(tq1^2-tq2^2)*c_3Psi ...
   +2*tq1*tq2*s_3Psi)+3*tau12tau22*((tq1*eps3-tq2*eps2)+4*(tq1*s_2Psi-tq2*c_2Psi) ...
   +(tq1^2-tq2^2)*s_3Psi-2*tq1*tq2*c_3Psi)+2*(3-2*eta_2)*(p1*tau2-p2*tau1) ...
   +2*(9-2*eta_2)*((p1^2-p2^2)*s_3Psi-2*p1*p2*c_3Psi) );
D_sp2_33 = ( tq2*(4-sig1_2)/(2*a2eta6sig22) )*( 4*tq1*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   +eta_2*(3*(p1*tau2-p2*tau1)+(p1^2-p2^2)*s_3Psi-2*p1*p2*c_3Psi) ) ...
   -( 3*tau12tau22/(2*a2eta6sig22) )*( 2*tq1*(10*tq1+3*(tq1*eps2+tq2*eps3) ...
   +6*(tq1*c_2Psi+tq2*s_2Psi)+(tq1^2-tq2^2)*c_3Psi+2*tq1*tq2*s_3Psi)+eta_2*(5+3*eps2+3*c_2Psi ...
   +(tq1*c_3Psi+tq2*s_3Psi)) ) - ( 2*tq1/(a2eta6sig22) )*( 3*(3-eta_2)*(p1*tau1+p2*tau2) ...
   +(9-eta_2)*((p1^2-p2^2)*c_3Psi+2*p1*p2*s_3Psi) );
D_sp2_34 = ( (4-sig1_2)/(2*a2eta6sig22) )*( (eta_2+4*tq2^2)*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   +eta_2*tq2*(3*(p1*tau1+p2*tau2)-(p1^2-p2^2)*c_3Psi-2*p1*p2*s_3Psi) ) ...
   - ( 3*tau12tau22/(2*a2eta6sig22) )*( 2*tq2*(10*tq1+3*(tq1*eps2+tq2*eps3) ...
   +6*(tq1*c_2Psi+tq2*s_2Psi)+(tq1^2-tq2^2)*c_3Psi+2*tq1*tq2*s_3Psi)+eta_2*(3*eps3+3*s_2Psi ...
   +(tq1*s_3Psi-tq2*c_3Psi)) ) - ( 2*tq2/(a2eta6sig22) )*( 3*(3-eta_2)*(p1*tau1+p2*tau2) ...
   +(9-eta_2)*((p1^2-p2^2)*c_3Psi+2*p1*p2*s_3Psi) );
D_sp2_35 = -( 1/(2*a2eta4sig23) )*( 2*tq2*p1*(9-sig1_2)*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   -2*tq2*sig2*(4-sig1_2)*(3*(xi1*s_Psi+xi2*c_Psi)+3*(p1*s_2Psi-p2*c_2Psi) ...
   +xj1*s_3Psi-xj2*c_3Psi) +3*(p1*c_2Psi*(sig3+4*p2^2) ...
   +p2*s_2Psi*(sig2-4*p1^2))*(10*tq1+3*(tq1*eps2+tq2*eps3)+6*(tq1*c_2Psi+tq2*s_2Psi) ...
   +(tq1^2-tq2^2)*c_3Psi+2*tq1*tq2*s_3Psi) +6*(3-2*eta_2)*(p1*c_Psi*(sig3+4*p2^2) ...
   +p2*s_Psi*(sig2-4*p1^2)) +2*(9-2*eta_2)*(p1*c_3Psi*(sig3+4*p2^2)+p2*s_3Psi*(sig2-4*p1^2)) );
D_sp2_36 = -( 1/(2*a2eta4sig23) )*( 2*tq2*p2*(9-sig1_2)*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   +2*tq2*sig2*(4-sig1_2)*(3*(xi1*c_Psi-xi2*s_Psi)+3*(p1*c_2Psi+p2*s_2Psi) ...
   +xj1*c_3Psi+xj2*s_3Psi) +3*(p1*s_2Psi*(sig2-4*p2^2) ...
   +p2*c_2Psi*(sig3+4*p1^2))*(10*tq1+3*(tq1*eps2+tq2*eps3)+6*(tq1*c_2Psi+tq2*s_2Psi) ...
   +(tq1^2-tq2^2)*c_3Psi+2*tq1*tq2*s_3Psi) +6*(3-2*eta_2)*(p1*s_Psi*(sig2-4*p2^2) ...
   -p2*c_Psi*(sig3+4*p1^2)) +2*(9-2*eta_2)*(p1*s_3Psi*(sig2-4*p2^2)-p2*c_3Psi*(sig3+4*p1^2)) );

D_sp2_41 = -(2/a) * tq2_sp2;
D_sp2_42 = -( 3/(4*a2eta4sig22) )*( 2*tq1*(4-sig1_2)*(2*tau12tau22+(tau1*xi1-tau2*xi2) ...
   +(p1^2-p2^2)*(tq1*c_3Psi+tq2*s_3Psi)+2*p1*p2*(tq1*s_3Psi-tq2*c_3Psi)) ...
   -4*tau1*tau2*(10*tq2-3*(tq1*eps3-tq2*eps2)+6*(tq1*s_2Psi-tq2*c_2Psi)+(tq1^2-tq2^2)*s_3Psi ...
   -2*tq1*tq2*c_3Psi)-3*tau12tau22*((tq1*eps2+tq2*eps3)-4*(tq1*c_2Psi+tq2*s_2Psi) ...
   -(tq1^2-tq2^2)*c_3Psi-2*tq1*tq2*s_3Psi)-2*(3-2*eta_2)*(p1*tau1+p2*tau2) ...
   +2*(9-2*eta_2)*((p1^2-p2^2)*c_3Psi+2*p1*p2*s_3Psi) );
D_sp2_43 = -( (4-sig1_2)/(2*a2eta6sig22) )*( (eta_2+4*tq1^2)*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   +eta_2*tq1*(3*(p1*tau2-p2*tau1)+(p1^2-p2^2)*s_3Psi-2*p1*p2*c_3Psi) ) ...
   - ( 3*tau12tau22/(2*a2eta6sig22) )*( 2*tq1*(10*tq2-3*(tq1*eps3-tq2*eps2) ...
   +6*(tq1*s_2Psi-tq2*c_2Psi)+(tq1^2-tq2^2)*s_3Psi-2*tq1*tq2*c_3Psi)-eta_2*(3*eps3-3*s_2Psi ...
   -(tq1*s_3Psi-tq2*c_3Psi)) ) + ( 2*tq1/(a2eta6sig22) )*( 3*(3-eta_2)*(p1*tau2-p2*tau1) ...
   -(9-eta_2)*((p1^2-p2^2)*s_3Psi-2*p1*p2*c_3Psi) );
D_sp2_44 = -( tq1*(4-sig1_2)/(2*a2eta6sig22) )*( 4*tq2*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   +eta_2*(3*(p1*tau1+p2*tau2)-(p1^2-p2^2)*c_3Psi-2*p1*p2*s_3Psi) ) ...
   - ( 3*tau12tau22/(2*a2eta6sig22) )*( 2*tq2*(10*tq2-3*(tq1*eps3-tq2*eps2) ...
   +6*(tq1*s_2Psi-tq2*c_2Psi)+(tq1^2-tq2^2)*s_3Psi-2*tq1*tq2*c_3Psi)+eta_2*(5+3*eps2-3*c_2Psi ...
   -(tq1*c_3Psi+tq2*s_3Psi)) ) + ( 2*tq2/(a2eta6sig22) )*( 3*(3-eta_2)*(p1*tau2-p2*tau1) ...
   -(9-eta_2)*((p1^2-p2^2)*s_3Psi-2*p1*p2*c_3Psi) );
D_sp2_45 = ( 1/(2*a2eta4sig23) )*( 2*tq1*p1*(9-sig1_2)*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   -2*tq1*sig2*(4-sig1_2)*(3*(xi1*s_Psi+xi2*c_Psi)+3*(p1*s_2Psi-p2*c_2Psi) ...
   +xj1*s_3Psi-xj2*c_3Psi) -3*(p1*c_2Psi*(sig3+4*p2^2) ...
   +p2*s_2Psi*(sig2-4*p1^2))*(10*tq2-3*(tq1*eps3-tq2*eps2)+6*(tq1*s_2Psi-tq2*c_2Psi) ...
   +(tq1^2-tq2^2)*s_3Psi-2*tq1*tq2*c_3Psi) +6*(3-2*eta_2)*(p1*s_Psi*(sig3+4*p2^2) ...
   -p2*c_Psi*(sig2-4*p1^2)) -2*(9-2*eta_2)*(p1*s_3Psi*(sig3+4*p2^2)-p2*c_3Psi*(sig2-4*p1^2)) );
D_sp2_46 = ( 1/(2*a2eta4sig23) )*( 2*tq1*p2*(9-sig1_2)*(6*tau1*tau2+3*(tau1*xi2+tau2*xi1) ...
   +(p1^2-p2^2)*(tq1*s_3Psi-tq2*c_3Psi)-2*p1*p2*(tq1*c_3Psi+tq2*s_3Psi)) ...
   +2*tq1*sig2*(4-sig1_2)*(3*(xi1*c_Psi-xi2*s_Psi)+3*(p1*c_2Psi+p2*s_2Psi) ...
   +xj1*c_3Psi+xj2*s_3Psi) -3*(p1*s_2Psi*(sig2-4*p2^2) ...
   +p2*c_2Psi*(sig3+4*p1^2))*(10*tq2-3*(tq1*eps3-tq2*eps2)+6*(tq1*s_2Psi-tq2*c_2Psi) ...
   +(tq1^2-tq2^2)*s_3Psi-2*tq1*tq2*c_3Psi) -6*(3-2*eta_2)*(p1*c_Psi*(sig2-4*p2^2) ...
   +p2*s_Psi*(sig3+4*p1^2)) +2*(9-2*eta_2)*(p1*c_3Psi*(sig2-4*p2^2)+p2*s_3Psi*(sig3+4*p1^2)) );

D_sp2_51 = -(2/a) * p1_sp2;
D_sp2_52 = ( 3*sig3*(1+eps2)/(2*a2eta4sig2) )*( (p1*s_2Psi-p2*c_2Psi) );
D_sp2_53 = -( sig3/(4*a2eta6sig2) )*( 4*tq1*(3*(xi1*c_Psi-xi2*s_Psi)+3*(p1*c_2Psi+p2*s_2Psi) ...
   +xj1*c_3Psi+xj2*s_3Psi)+eta_2*(3*tau1+(p1*c_3Psi+p2*s_3Psi)) );
D_sp2_54 = -( sig3/(4*a2eta6sig2) )*( 4*tq2*(3*(xi1*c_Psi-xi2*s_Psi)+3*(p1*c_2Psi+p2*s_2Psi) ...
   +xj1*c_3Psi+xj2*s_3Psi)-eta_2*(3*tau2-(p1*s_3Psi-p2*c_3Psi)) );
D_sp2_55 = ( 1/(4*a2eta4sig22) )*( 4*p1*(3*(xi1*c_Psi-xi2*s_Psi)+3*(p1*c_2Psi+p2*s_2Psi) ...
   +xj1*c_3Psi+xj2*s_3Psi)-sig2*sig3*(3*(tq1*c_Psi-tq2*s_Psi) ...
   +3*c_2Psi+(tq1*c_3Psi+tq2*s_3Psi)) );
D_sp2_56 = ( 1/(4*a2eta4sig22) )*( 4*p2*(3*(xi1*c_Psi-xi2*s_Psi)+3*(p1*c_2Psi+p2*s_2Psi) ...
   +xj1*c_3Psi+xj2*s_3Psi)-sig2*sig3*(3*(tq1*s_Psi+tq2*c_Psi) ...
   +3*s_2Psi+(tq1*s_3Psi-tq2*c_3Psi)) );

D_sp2_61 = -(2/a) * p2_sp2;
D_sp2_62 = -( 3*sig3*(1+eps2)/(2*a2eta4sig2) )*( (p1*c_2Psi+p2*s_2Psi) );
D_sp2_63 = -( sig3/(4*a2eta6sig2) )*( 4*tq1*(3*(xi1*s_Psi+xi2*c_Psi)+3*(p1*s_2Psi-p2*c_2Psi) ...
   +xj1*s_3Psi-xj2*c_3Psi)+eta_2*(3*tau2+(p1*s_3Psi-p2*c_3Psi)) );
D_sp2_64 = -( sig3/(4*a2eta6sig2) )*( 4*tq2*(3*(xi1*s_Psi+xi2*c_Psi)+3*(p1*s_2Psi-p2*c_2Psi) ...
   +xj1*s_3Psi-xj2*c_3Psi)+eta_2*(3*tau1-(p1*c_3Psi+p2*s_3Psi)) );
D_sp2_65 = ( 1/(4*a2eta4sig22) )*( 4*p1*(3*(xi1*s_Psi+xi2*c_Psi)+3*(p1*s_2Psi-p2*c_2Psi) ...
   +xj1*s_3Psi-xj2*c_3Psi)-sig2*sig3*(3*(tq1*s_Psi+tq2*c_Psi) ...
   +3*s_2Psi+(tq1*s_3Psi-tq2*c_3Psi)) );
D_sp2_66 = ( 1/(4*a2eta4sig22) )*( 4*p2*(3*(xi1*s_Psi+xi2*c_Psi)+3*(p1*s_2Psi-p2*c_2Psi) ...
   +xj1*s_3Psi-xj2*c_3Psi)+sig2*sig3*(3*(tq1*c_Psi-tq2*s_Psi) ...
   +3*c_2Psi+(tq1*c_3Psi+tq2*s_3Psi)) );

D_sp2 = [ D_sp2_11  D_sp2_12  D_sp2_13  D_sp2_14  D_sp2_15  D_sp2_16;
          D_sp2_21  D_sp2_22  D_sp2_23  D_sp2_24  D_sp2_25  D_sp2_26;
          D_sp2_31  D_sp2_32  D_sp2_33  D_sp2_34  D_sp2_35  D_sp2_36;
          D_sp2_41  D_sp2_42  D_sp2_43  D_sp2_44  D_sp2_45  D_sp2_46;
          D_sp2_51  D_sp2_52  D_sp2_53  D_sp2_54  D_sp2_55  D_sp2_56;
          D_sp2_61  D_sp2_62  D_sp2_63  D_sp2_64  D_sp2_65  D_sp2_66 ];

% Evaluating equinoctial variables
a_osc = a + coef * ( a_lp + a_sp1 + a_sp2 );
Psi_osc = Psi + coef * ( Psi_lp + Psi_sp1 + Psi_sp2 );

tq1_osc = tq1 + coef * ( tq1_lp + tq1_sp1 + tq1_sp2 );
tq2_osc = tq2 + coef * ( tq2_lp + tq2_sp1 + tq2_sp2 );
p1_osc = p1 + coef * ( p1_lp + p1_sp1 + p1_sp2 );
p2_osc = p2 + coef * ( p2_lp + p2_sp1 + p2_sp2 );

% Transformation Matrix D_J2
D_J2 = DI + coef * ( D_lp + D_sp1 + D_sp2 );

% Osculating variables from Mean variables
osc_equi = [ a_osc;  Psi_osc;  tq1_osc;  tq2_osc;  p1_osc;  p2_osc ];

end