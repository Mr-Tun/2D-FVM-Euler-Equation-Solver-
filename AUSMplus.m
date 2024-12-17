function [F] = AUSMplus(UL,UR,n)
%AUSM function is used to compute the flux by AUSM method
gamma = 1.4;
%% decouple the normal vector
nx = n(1);
ny = n(2);
%% decouple the solution vector U left
rhol = UL(1);
rhoul = UL(2);
rhovl = UL(3);
rhoel = UL(4);
ul = rhoul / rhol;
vl = rhovl / rhol;
unl = nx*ul + vl*ny;
utl = -ul*ny + vl*nx;
el = rhoel / rhol - 0.5 * (ul^2+vl^2);
pl = (gamma - 1) * rhol * el;
al =  (gamma.*pl./rhol)^0.5;
Hl = (rhoel+pl)/rhol;
%% decouple the solution vector U right
rhor = UR(1);
rhour = UR(2);
rhovr = UR(3);
rhoer = UR(4);
ur = rhour / rhor;
vr = rhovr / rhor;
unr = nx*ur + ny*vr;
utr = -ur*ny + vr*nx;
er = rhoer ./ rhor - 0.5*(ur^2+vr^2);
pr = (gamma - 1) * rhor* er;
ar =  (gamma.*pr./rhor)^0.5;
Hr = (rhoer+pr)/rhor;
%% AUSM Method which support 1st and 2nd order reconstruction
Mal = unl/al;
Mar = unr/ar;
Ma_plus = Mal*(Mal>1)+0.25*(Mal+1)^2*(abs(Mal)<=1);
Ma_minus = Mar*(Mar<-1)-0.25*(Mar-1)^2*(abs(Mar)<=1);
p_plus = pl*(Mal>1)+pl*(1+Mal)^2/4*(2-Mal)*(abs(Mal)<=1);
p_minus = pr*(Mar<-1)+pr*(1-Mar)^2/4*(2+Mar)*(abs(Mar)<=1);
p_mid = p_plus + p_minus;
fc_plus = rhol*al*Ma_plus*[1, ul,vl,Hl];
fc_minus = rhor*ar*Ma_minus*[1,ur,vr,Hr];
Fc_mid = fc_minus +fc_plus;
Fc_mid(2) = Fc_mid(2)+p_mid*nx;
Fc_mid(3) = Fc_mid(3)+p_mid*ny;
F = Fc_mid;
%% AUSM+ Method, which only support 1st order reconstruction
% a_mid = (al+ar)/2;
% Mal = unl/a_mid;
% Mar = unr/a_mid;
% M_mid = Mal*(Mal>1)+0.25*(Mal+1)^2*(abs(Mal)<=1) + Mar*(Mar<-1)-0.25*(Mar-1)^2*(abs(Mar)<=1);
% u_mid = M_mid*a_mid;
% phil = UL; phir = UR;
% phil(4) = Hl*rhol; phir(4) = rhor*Hr;
% phi_mid = phil*(u_mid>=0)+phir*(u_mid<0);
% Fc_mid = phi_mid * u_mid;
% p_mid = pl*(1*(Mal>1)+(2-Mal)*0.25*(Mal+1)^2*(abs(Mal)<=1)) + pr*(1*(Mar<-1)+(2+Mar)*0.25*(Mar-1)^2*(abs(Mar)<=1));
% Fc_mid(2) = Fc_mid(2)+p_mid*nx;
% Fc_mid(3) = Fc_mid(3)+p_mid*ny;
% F = Fc_mid;
end
 