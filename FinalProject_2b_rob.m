% Script: convect1d.m
clear all;

%set constants
v_max = 33;     % m/s
p_max = 0.25;   % 1/m
p_L = 0.15;     % Left boundary density
p_R = 0.25;     % Right boundary density

% Set-up grid
xL = 0;
xR = 2000;
Nx = 100;       % number of cells
x = linspace(xL,xR,Nx+1);

% Calculate midpoint of each cell. 
xmid = zeros(1,Nx);
for i = 1:1:Nx
    xmid(i) = (.5)*(x(i) + x(i+1));
end % for i = 1:1:Nx

% Calculate cell size. We assume that each cell is the same size
dx = x(2) - x(1);

% Set shock velocity
a = ((1-p_L/p_max)*v_max*p_L - (1-p_R/p_max)*v_max*p_R)/(p_L-p_R);

% Set final time
tfinal = 35;

% Set timestep
CFL = 0.5;
dt = CFL*dx/v_max;
t=0;

% Set initial condition to p0 = 0.15 (1/m)
P = zeros(1,numel(xmid));       % Density in each cell
for i = 1:1:numel(xmid)
    P(i) = .15;
end % for i = 1:1:numel(xmid)

% Initialize exact solution array
P_Exact=zeros(size(P));         % Exact density in each cell (at midpoints)

% Loop until t > tfinal
while (t < tfinal)
    Pbc = [p_L, P, p_R]; % This enforces the bc
    
    % Use P_j+(1/2) to caculate q_j+(1/2) at each boundary 
    Q_j_half = zeros(1,Nx+2);
    for i = 2:1:(Nx+1)
        s1 = (1-2*Pbc(i)/p_max)*v_max;
        s2 = (1-2*Pbc(i+1)/p_max)*v_max;
        if s1>=s2
            s_max=s1;
        else
            s_max=s2;
        end
        
        % Note: Pbc(i) = P(i-1) for i in [2,Nx+1]
        Q_j_half(i-1) = 0.5*(Pbc(i)*(1-Pbc(i)/p_max)*v_max+Pbc(i+1)...
            *(1-Pbc(i+1)/p_max)*v_max)+ (s_max/2)*(Pbc(i)-Pbc(i+1));
    end % for i = 1:1:(Nx+1)
    
    % Update P for all interior elements of Pbc
    for i = 1:1:Nx
        P(i) = P(i) + (dt/dx)*(Q_j_half(i) - Q_j_half(i+1));
    end % for i = 2:1:(Nx+2)
    
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
    figure(1);
    clf;    
    hold on;
    grid on;

    plot(x,[P, P(Nx)],'go');
    plot(x,[P_Exact, P_Exact(Nx)],'k-');
        axis([xL, xR, 0, 0.3]);
        legend('FVM','Exact');
        title(['t=',num2str(t)]);
        ylabel('density (1/m)');
        xlabel('distance (m)');
    drawnow;
    hold off;
end % while(t < tfinal)