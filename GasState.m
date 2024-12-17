function [rho] = GasState(p,T)
%GasState function uses the ideal gas equation of state to solve the air's temperature by
%pressure and density 
% global Rg
Rg = 287.05;
rho = p/Rg/T;  % temperature of the gas --k
end
