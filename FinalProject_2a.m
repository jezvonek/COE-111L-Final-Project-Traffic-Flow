% Script: convect1d.m
clear all;

%set constants
v_max = 33; %m/s
p_max = 0.25; %1/m
p_L = 0.15;
p_R = 0.25;

% Set-up grid
xL = 0;
xR = 2000;
Nx = 100; % number of control volumes
x = linspace(xL,xR,Nx+1);
% Calculate midpoint values of x in each control volume
xmid = 0.5*(x(1:Nx) + x(2:Nx+1));

% Calculate cell size in control volumes (assumed equal)
dx = x(2) - x(1);
% Set shock velocity
a = ((1-p_L/p_max)*v_max*p_L - (1-p_R/p_max)*v_max*p_R)/(p_L-p_R);
% Set final time
tfinal = 35;
% Set timestep
CFL = 0.5;
dt = CFL*dx/v_max;
t=0;

% Set initial condition to Q0 = 0.15 1/m
P = 0.15*ones(size(xmid)); %1/m

%initialize size of exact solution array
P_Exact=zeros(size(P));

% Loop until t > tfinal
while (t < tfinal)
Pbc = [p_L, P, p_R]; % This enforces the bc
% Calculate the flux at each interface
%Pbc_avg = 0.5*(Pbc(1:Nx+1)+Pbc(2:Nx+2));
Q = (1-Pbc/p_max)*v_max.*Pbc;
% Calculate rhs in each cell
R = Q(3:Nx+2) - Q(2:Nx+1);
% Forward Euler step
P = P - (dt/dx)*R;

%calculate exact solution at this time step
    for i=1:length(P_Exact)
        if ((dx*(i-1/2))-2000)/t < a
            P_Exact(i) = p_L;
        else
            P_Exact(i) = p_R;
        end
    end
    
% Increment time
t = t + dt;

% Plot current solution
figure(1)
clf
plot(x,[P, P(Nx)],'go');
hold on;
plot(x,[P_Exact, P_Exact(Nx)],'k-');
axis([xL, xR, 0, 0.3]);
legend('FVM','Exact');
title(['t=',num2str(t)]);
ylabel('density (1/m)');
xlabel('distance (m)');
grid on;
drawnow;
hold off;
end
