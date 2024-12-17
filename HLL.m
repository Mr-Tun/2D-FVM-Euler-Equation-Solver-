function [F] = HLL(UL,UR,n)
% HLL function is used to get the approximate Riemann solver by the HLL
% scheme and make 2D Riemann problem to 1D by transferring coorindate
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
%% compute Roe averge variable to get more precise wave speed
RT = sqrt(rhor/rhol);
u_Roe = (ul+RT*ur)/(1+RT);
v_Roe = (vl+RT*vr)/(1+RT);
H_Roe = (Hl+RT*Hr)/(1+RT);
a_Roe = sqrt((gamma-1)*(H_Roe - (u_Roe*u_Roe + v_Roe*v_Roe)/2));
un_Roe = u_Roe*nx + v_Roe*ny;
%% HLL Scheme 
FL = [rhol*unl,rhoul*unl+pl*nx,rhovl*unl+pl*ny,rhoel*unl+pl*unl];
FR = [rhor*unr,rhour*unr+pr*nx,rhovr*unr+pr*ny,rhoer*unr+pr*unr];
%     FL = [rhol*unl,rhol*unl*unl+pl,rhol*utl*unl,rhoel*unl+pl*unl];
%     FR = [rhor*unr,rhor*unr*unr+pr,rhor*utr*unr,rhoer*unr+pr*unr];
% UL(2) = rhol*unl;    UL(3) = rhol*utl;
% UR(2) = rhor*unr;    UR(3) = rhor*utr;
 SL = min([unl-al,unr-ar]);
 SR = max([unl+al,unr+ar]);

     % SL = min([unl-al,un_Roe-a_Roe]);
     % SR = max([unl+al,un_Roe+a_Roe]);
if SL>=0
    F = FL;
elseif SR<=0
    F = FR;
else
    F = (SR*FL -SL*FR + SL*SR*(UR-UL))./(SR-SL);
end
%recover the flux
%       F(2) = F(2)*nx-F(3)*ny;
%       F(3) = F(2)*ny+F(3)*nx;
end

