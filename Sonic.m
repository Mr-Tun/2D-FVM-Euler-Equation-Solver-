function [a] = Sonic(T)
%Sonic function is used to compute the local sonic velocity
gamma = 1.4;
Rg = 287.05;
a = sqrt(gamma*Rg*T);
end

