function [rho, u, v, T] = Boundary(rho, u, v, T, lower)
% Boundary function is used to set the boundary condition on left, right
% and lower condition by using ghost cells

%% Lelf & right
rho(2,3:end-2) = rho(3,3:end-2); rho(end-1,3:end-2)=rho(end-2,3:end-2);
rho(1,3:end-2) = rho(3,3:end-2); rho(end,3:end-2)=rho(end-2,3:end-2);
u(2,3:end-2) = u(3,3:end-2); u(end-1,3:end-2)=u(end-2,3:end-2);
u(1,3:end-2) = u(3,3:end-2); u(end,3:end-2)=u(end-2,3:end-2);
v(2,3:end-2) = v(3,3:end-2); v(end-1,3:end-2)=v(end-2,3:end-2);
v(1,3:end-2) = v(3,3:end-2); v(end,3:end-2)=v(end-2,3:end-2);
T(2,3:end-2) = T(3,3:end-2); T(end-1,3:end-2)=T(end-2,3:end-2);
T(1,3:end-2) = T(3,3:end-2); T(end,3:end-2)=T(end-2,3:end-2);
%% lower
lowerx = squeeze(lower(:,:,1));
lowery = squeeze(lower(:,:,2));

rho(3:end-2,2) = rho(3:end-2,3); rho(3:end-2,1) = rho(3:end-2,4); 

un_ghost1 = -(u(3:end-2,3).*lowerx(:,1)+v(3:end-2,3).*lowery(:,1)); 
un_ghost2 = -(u(3:end-2,4).*lowerx(:,2)+v(3:end-2,4).*lowery(:,2)); 
ut_ghost1 = -u(3:end-2,3).*lowery(:,1)+v(3:end-2,3).*lowerx(:,1); 
ut_ghost2 = -u(3:end-2,4).*lowery(:,2)+v(3:end-2,4).*lowerx(:,2);

T(3:end-2,2) = T(3:end-2,3); T(3:end-2,1) = T(3:end-2,4); 

u(3:end-2,2) = un_ghost1.*lowerx(:,1) - ut_ghost1.*lowery(:,1);
u(3:end-2,1) = un_ghost2.*lowerx(:,2) - ut_ghost2.*lowery(:,2);
v(3:end-2,2) = un_ghost1.*lowery(:,1) + ut_ghost1.*lowerx(:,1);
v(3:end-2,1) = un_ghost2.*lowery(:,2) + ut_ghost2.*lowerx(:,2);

end

