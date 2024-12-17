clc; clear;close all;
%% Define the parameter of the test gas (air)
% gobal gamma Rg Cv
gamma=1.4;                      % the adiabatic index γ for air is approximately 1.4
Rg=287.05;                      % the specific gas constant for air --J/(kg·K)
Cv = 1/(gamma-1)*Rg;            % compute the specific heat capacity at constant volume  --J/(kg·K)
%% Read the grid of the case
[X,Y]=ReadPLOT3D('Cylinder.dat');
[nx, ny] = size(X);
[Xc,Yc,omega,SI,SJ,right,upper,left,lower] = MeshSetting(X,Y);
%% define the flow conditions and initialized the grid
M0 = 8.1;                         % the Ma number of inflow
T0 = 63.73;                       % the temperature of inflow  -K
p0 = 370.6;                       % the pressure of inflow -Pa
rho0 = GasState(p0,T0);           % compute the density of inflow -kg/m³
u0 = M0 * Sonic(T0);              % compute the velocity of the inflow -m/s
v0 = 0;                           % initialize the velocity in y-axis
[p_total] = NewtonianLaw(M0,p0,rho0,u0,Xc,Yc,lower); % predict the maximum pressure
EmpiricalShockShape(0.02,M0);     % predict the empirical shock shape
rho = rho0*ones(nx+3, ny+3);          % initial the primitive variables
rho(2:end-2,end-1:end) = rho0; rho(1:2,1:2) = NaN; rho(1:2,end-1:end) = NaN; rho(end-1:end,end-1:end) =NaN; rho(end-1:end,1:2) =NaN;
u = u0*ones(nx+3,ny+3); u(2:end-2,end-1:end) = u0; u(1:2,1:2) = NaN; u(1:2,end-1:end) = NaN; u(end-1:end,end-1:end) =NaN; u(end-1:end,1:2) =NaN;
v = v0*ones(nx+3,ny+3); v(2:end-2,end-1:end) = v0; v(1:2,1:2) = NaN; v(1:2,end-1:end) = NaN; v(end-1:end,end-1:end) =NaN; v(end-1:end,1:2) =NaN;
p = p0*ones(nx+3,ny+3); p(2:end-2,end-1:end) = p0; p(1:2,1:2) = NaN; p(1:2,end-1:end) = NaN; p(end-1:end,end-1:end) =NaN; p(end-1:end,1:2) =NaN;
T = T0*ones(nx+3,ny+3); T(2:end-2,end-1:end) = T0; T(1:2,1:2) = NaN; T(1:2,end-1:end) = NaN; T(end-1:end,end-1:end) =NaN; T(end-1:end,1:2) =NaN;
rho_new = rho; u_new=u; v_new=v; T_new=T;
t = 0;
timestep = 0;                      % initial timestep
time = [];                         % initial history vector of time
time(1,1)=t;
residual = 1.1;
inflow_rho = rho0*ones(nx-1,2); inflow_u = u0*ones(nx-1,2);
inflow_v = v0*ones(nx-1,2); inflow_p = p0*ones(nx-1,2);
inflow_T = T0*ones(nx-1,2);
min_sj=min(SJ(:)); min_si=min(SI(:));
min_dx=min(min_si,min_sj);
Fleft = zeros(nx-1,ny-1,4); Fright = zeros(nx-1,ny-1,4); Fupper = zeros(nx-1,ny-1,4); Flower = zeros(nx-1,ny-1,4);
%% Main Loop
while residual>0.05&&t<0.0005                        % set end conditions                                 % initial time --s
    [rho, u, v, T] = Boundary(rho, u, v, T, lower); % set the boundary condition before time marching
    % velocity_max = abs(u) + (gamma*Rg*T).^0.5;
    % dt = 0.1*min_dx/max(velocity_max(:));           % CFL condition
    dt = 1*10^(-8);                                 %CFL number roughly to 0.2
    t = t+dt;                                       % current time in shock tube
    timestep = timestep+1;                          % current timestep
    if timestep ==1
        disp('the FVM solver is begin to compute, pls wait for a while');
    end
    time(timestep)=t;
    for i=3:nx+1
        for j=3:ny+1
            % transfer the primitive variables to solution vector U
            U_I_J = [rho(i,j),rho(i,j)*u(i,j),rho(i,j)*v(i,j),rho(i,j)*(Cv*T(i,j)+0.5*u(i,j)^2+0.5*v(i,j)^2)];
            U_I_1_J =[rho(i-1,j),rho(i-1,j)*u(i-1,j),rho(i-1,j)*v(i-1,j),rho(i-1,j)*(Cv*T(i-1,j)+0.5*u(i-1,j)^2+0.5*v(i-1,j)^2)];
            U_I_2_J =[rho(i-2,j),rho(i-2,j)*u(i-2,j),rho(i-2,j)*v(i-2,j),rho(i-2,j)*(Cv*T(i-2,j)+0.5*u(i-2,j)^2+0.5*v(i-2,j)^2)];
            U_I1_J = [rho(i+1,j),rho(i+1,j)*u(i+1,j),rho(i+1,j)*v(i+1,j),rho(i+1,j)*(Cv*T(i+1,j)+0.5*u(i+1,j)^2+0.5*v(i+1,j)^2)];
            U_I2_J = [rho(i+2,j),rho(i+2,j)*u(i+2,j),rho(i+2,j)*v(i+2,j),rho(i+2,j)*(Cv*T(i+2,j)+0.5*u(i+2,j)^2+0.5*v(i+2,j)^2)];
            U_I_J_1 =[rho(i,j-1),rho(i,j-1)*u(i,j-1),rho(i,j-1)*v(i,j-1),rho(i,j-1)*(Cv*T(i,j-1)+0.5*u(i,j-1)^2+0.5*v(i,j-1)^2)];
            U_I_J_2 =[rho(i,j-2),rho(i,j-2)*u(i,j-2),rho(i,j-2)*v(i,j-2),rho(i,j-2)*(Cv*T(i,j-2)+0.5*u(i,j-2)^2+0.5*v(i,j-2)^2)];
            U_I_J1 = [rho(i,j+1),rho(i,j+1)*u(i,j+1),rho(i,j+1)*v(i,j+1),rho(i,j+1)*(Cv*T(i,j+1)+0.5*u(i,j+1)^2+0.5*v(i,j+1)^2)];
            U_I_J2 = [rho(i,j+2),rho(i,j+2)*u(i,j+2),rho(i,j+2)*v(i,j+2),rho(i,j+2)*(Cv*T(i,j+2)+0.5*u(i,j+2)^2+0.5*v(i,j+2)^2)];
            % compute F_left at i-2,j-2 cell
            [Uleft,Uright] = Reconstruction(U_I1_J,U_I_J,U_I_1_J, U_I_2_J);
            Fleft(i-2,j-2,:)= HLL(Uleft, Uright, left(i-2,j-2,:));
            % Fleft(i-2,j-2,:)= AUSMplus(Uleft, Uright, left(i-2,j-2,:));
            % compute F_right at i-2,j-2 cell
            [Uleft,Uright] = Reconstruction(U_I_1_J,U_I_J,U_I1_J,U_I2_J);
            Fright(i-2,j-2,:) = HLL(Uleft, Uright, right(i-2,j-2,:));
            % Fright(i-2,j-2,:) = AUSMplus(Uleft, Uright, right(i-2,j-2,:));
            % compute F_lower at i-2,j-2 cell
            [Ulower,Uupper] = Reconstruction(U_I_J1,U_I_J,U_I_J_1,U_I_J_2);
            Flower(i-2,j-2,:) = HLL(Ulower,Uupper, lower(i-2,j-2,:));
            % Flower(i-2,j-2,:) = AUSMplus(Ulower,Uupper, lower(i-2,j-2,:));
            % compute F_upper at i-2,j-2 cell
            [Ulower,Uupper] = Reconstruction(U_I_J_1,U_I_J,U_I_J1,U_I_J2);
            Fupper(i-2,j-2,:) = HLL(Ulower,Uupper, upper(i-2,j-2,:));
            % Fupper(i-2,j-2,:) = AUSMplus(Ulower,Uupper, upper(i-2,j-2,:));
        end
    end
    for i=3:nx+1
        for j=3:ny+1
            % compute the Flux sum at i-2,j-2 cell and update U at next timestep
            F_left = [Fleft(i-2,j-2,1),Fleft(i-2,j-2,2),Fleft(i-2,j-2,3),Fleft(i-2,j-2,4)];
            F_right = [Fright(i-2,j-2,1),Fright(i-2,j-2,2),Fright(i-2,j-2,3),Fright(i-2,j-2,4)];
            F_lower = [Flower(i-2,j-2,1),Flower(i-2,j-2,2),Flower(i-2,j-2,3),Flower(i-2,j-2,4)];
            F_upper = [Fupper(i-2,j-2,1),Fupper(i-2,j-2,2),Fupper(i-2,j-2,3),Fupper(i-2,j-2,4)];
            sum = F_left*SI(i-2,j-2)+F_right*SI(i+1-2,j-2)+F_lower*SJ(i-2,j-2)+F_upper*SJ(i-2,j+1-2);
            sum = -sum/omega(i-2,j-2);
            U_next_time = dt * sum + [rho(i,j),rho(i,j)*u(i,j),rho(i,j)*v(i,j),rho(i,j)*(Cv*T(i,j)+0.5*u(i,j)^2+0.5*v(i,j)^2)];
            % decouple the U at i-2,j-2 and update the primitive variables
            [rho_new(i,j),u_new(i,j),v_new(i,j),T_new(i,j)] = Decouple(U_next_time);
        end
    end
    p(3:end-2,3:end-2)=rho(3:end-2,3:end-2).*Rg.*T(3:end-2,3:end-2);
    rho(3:end-2,3:end-2)=rho_new(3:end-2,3:end-2);
    u(3:end-2,3:end-2)=u_new(3:end-2,3:end-2);
    v(3:end-2,3:end-2)=v_new(3:end-2,3:end-2);
    T(3:end-2,3:end-2)=T_new(3:end-2,3:end-2);
    if mod(timestep, 500) == 0
        p_new = rho.*Rg.*T;
        error_p = abs(p_new(3:end-2,3:end-2)-p(3:end-2,3:end-2));
        residual = 0;
        for i = 1:numel(error_p)
            residual = residual + error_p(i);
        end
        residual = residual/(nx-1)/(ny-1);
        real_p = p_new(3:end-2, 3:end-2);
        p_max = max(real_p(:));
        fprintf('Timestep: %d, residual: %.8e, Max Pressure: %d, Current Time: %.8e\n', timestep, residual, p_max, t);
        figure(3)
        pcolor(Xc, Yc, p_new(3:end-2,3:end-2));
        shading interp;
        colorbar;
        cb = colorbar;
        cb.Ruler.Exponent = 0;
        colormap(jet);
        title('Pressure Distribution');
        xlabel('X/m');
        ylabel('Y/m');
        axis equal;
        hold on
    end
end
p = rho.*Rg.*T;
figure(2);hold on; plot(Yc(:,1),p(3:end-2,1));legend('Numerical Solution','Newtonian Law');
