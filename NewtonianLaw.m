function [p_total] = NewtonianLaw(Ma,p,rho,v,Xc,Yc,lower)
%NewtonianLaw is used to compute the estimated value of pressure
nx_wall = squeeze(lower(:,1,1)); ny_wall = squeeze(lower(:,1,2));
sita = atan(nx_wall./ny_wall);
Xc_wall = Xc(:,1); Yc_wall = Yc(:,1);
gamma = 1.4;
Cpmax = 2/(gamma*Ma^2)*(((gamma+1)^2*Ma^2/(4*gamma*Ma^2-2*(gamma-1)))^(gamma/(gamma-1))*((1-gamma+2*gamma*Ma^2)/(gamma+1))-1);
p_total=Cpmax*(0.5*rho*v*v)+p;
Cp = Cpmax*(sin(sita)).^2;
p_wall =0.5*rho*v^2*Cp+p;
figure(2)
plot(Yc_wall,p_wall);
legend('Newtonian Law')
xlabel('Y/m');
ylabel('Pressure/Pa');
title('The Pressure Distribution at Cylinder Wall')
ax = gca;
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
disp(['the total pressure at the center point of the cylinder wall is ',num2str(p_total),'Pa']);
end

