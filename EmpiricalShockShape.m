function [] = EmpiricalShockShape(R,Ma)
%EmpiricalShockShape is used to compute the shape of the shock 

sigma = 0.386*exp(4.67/Ma/Ma)*R;                 % compute the shock detachment distance
Rc = 1.386*exp(1.8/(Ma-1)^0.75)*R;               % compute the  radius of curvature of the shock wave at the vertex of the  hyperbola
beta = asin(1/Ma);                               % compute the angle of mach wave
disp(['Shock detachment distance σ = ',num2str(sigma*1000),'mm, Rc =',num2str(Rc),'mm, the angle of Mach wave β =', num2str(beta)]);

y = -2.5*R:2.5*R/160:2.5*R;                               
x = R+sigma-Rc*(1/tan(beta))^2* ((1+y.^2.*tan(beta).*tan(beta)./Rc./Rc).^0.5-1);% Billig relationships
x = -x;
figure(1)
hold on 
plot(x, y, 'LineWidth', 2,'color','black','LineStyle',':');
hold off
end

