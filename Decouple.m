function [rho,u,v,T] = Decouple(U)
%Decouple function is used to decouple the U to primitive variables

gamma=1.4;                      % the adiabatic index γ for air is approximately 1.4
Rg=287.05;                      % the specific gas constant for air --J/(kg·K)
Cv = 1/(gamma-1)*Rg;            % compute the specific heat capacity at constant volume  --J/(kg·K)
rho = U(1);                 %decouple the solution vector U
rhou = U(2);
rhov =U(3);
rhoe = U(4);
u = rhou / rho;                  % compute lambda matrix
v = rhov / rho;
e = rhoe / rho - 0.5 * (u ^ 2+v ^ 2);
T = e/Cv;

end

