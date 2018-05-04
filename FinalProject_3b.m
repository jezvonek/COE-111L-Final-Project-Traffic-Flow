% Script: convect1d.m
clear all;

%set constants
v_max = 33;     % m/s
p_max = 0.25;   % car/m
flux_l = 1.2;   % car/s
p_L = p_max/2 - .5*(p_max^2-4*p_max*flux_l/v_max)^(1/2);     % Left boundary density
p_R = 0.02;     % Right boundary density

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
tfinal = 20;
% Set timestep
CFL = 0.5;
dt = CFL*dx/v_max;
t=0;

% Set initial condition to p0 = 0.15 (1/m)
P = zeros(1,Nx);       % Density in each cell
for i = 1:1:Nx
    P(i) = .02;
end % for i = 1:1:numel(xmid)

%initialize size of exact solution array
P_Exact=zeros(size(P));
Q=zeros(1,length(P)+2);
% Loop until t > tfinal
while (t < tfinal)
Pbc = [p_L, P, p_R]; % This enforces the bc
% Calculate the flux at each interface
    for i=1:Nx+1
        s1 = (1-2*Pbc(i)/p_max)*v_max;
        s2 = (1-2*Pbc(i+1)/p_max)*v_max;
        if s1>=s2
            s_max=s1;
        else
            s_max=s2;
        end
        
        Q(i) = 0.5*(Pbc(i)*(1-Pbc(i)/p_max)*v_max+Pbc(i+1)*(1-Pbc(i+1)/p_max)*v_max)+ (s_max/2)*(Pbc(i)-Pbc(i+1));
    end
    
% Calculate rhs in each cell
R = Q(2:Nx+1) - Q(1:Nx);
% Forward Euler step
P = P - (dt/dx)*R;

%  %calculate exact solution at this time step
        for i=1:length(P_Exact)
          if (i-1/2)*dx/t <= (1-2*p_L/p_max)*v_max
              P_Exact(i) = p_L;
          elseif (1-2*p_L/p_max)*v_max < (i-1/2)*dx/t && (i-1/2)*dx/t < (1-2*p_R/p_max)*v_max
              P_Exact(i) = 0.5*p_max*(1-(i-1/2)*dx/(v_max*t));
          else
              P_Exact(i) = p_R;
          end
        end
    
% Increment time
t = t + dt;

% Plot current solution
figure(1)
clf
plot(x,[P, P(Nx)],'mo');
hold on;
plot(x,[P_Exact, P_Exact(Nx)],'k-');
axis([xL, xR, 0, 0.06]);
legend('FVM','Exact');
title(['t=',num2str(t)]);
ylabel('density (1/m)');
xlabel('distance (m)');
grid on;
drawnow;
hold off;
end