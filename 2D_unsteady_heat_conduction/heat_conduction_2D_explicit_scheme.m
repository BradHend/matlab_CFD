% Solution of 2-D Unsteady Heat Conduction equation using Explicit Euler
% scheme for time
%   PDE:
%   dT_dt = alpha*(d2T_dx2 + d2T_dy2)

clear, clc, close all

%grid sizes to evaluate (NxN)
grid_sizes = [10, 20, 40];
%% Set up input params. for solver routine
%domain params.
domain.L = 1.0;
%boundary conditions
boundary_conditions.T_x0y = 0.0; %left
boundary_conditions.dTdx_xLy = 0.0; %right
boundary_conditions.T_xy0 = 0.0; %bottom
boundary_conditions.T_xyL = 100.0; %top
boundary_conditions.T_0internal = 0.0; %internal points
%solution control
solver.dt = 1e-04;
solver.t_stop = 0.05;
% solver.t_stop = 50.0; %"steady-state"
solver.n_iterations = ceil(solver.t_stop/solver.dt); %figure out how many time iterations to do

%% calling Explicit Euler function, plotting results
figure(1)
figure(2)
for k = 1:numel(grid_sizes) %loop over all grids
    n_nodes = grid_sizes(k);
    %make grid coords. for plotting
    x = linspace(0,domain.L,n_nodes);
    y = linspace(0,domain.L,n_nodes);
    %call explicit euler solver for current n_nodes, etc.
    T = explicit_euler(n_nodes, domain, boundary_conditions, solver);
    %slice idxs for plotting
    dx = domain.L/n_nodes;
    dy = domain.L/n_nodes;
    x_plt_idx = 0.5*domain.L/dx;
    y_plt_idx = 0.5*domain.L/dy;
    figure(1) %plot T(x=0.5,y)
    hold on
    plot(y, T(x_plt_idx,:),'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])
    figure(2) %plot T(x,y=0.5)
    hold on
    plot(x, T(:,y_plt_idx),'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])
end %grid loop

figure(1) %Temperature slice @(x=0.5,y)
title(['T(x=0.5,y) @ t=',num2str(solver.t_stop),'s', ', dt=',num2str(solver.dt)])
legend('location','Best')
box on, grid on
ylabel('T')
xlabel('y coordinate')

figure(2) %Temperature slice @(x,y=0.5)
title(['T(x,y=0.5) @ t=',num2str(solver.t_stop),'s', ', dt=',num2str(solver.dt)])
legend('location','Best')
box on, grid on
ylabel('T')
xlabel('x coordinate')

figure(3) %temperature contour
contourf(x, y, T')
title(['Temperature Contour @t=',num2str(solver.t_stop),'s',': ',[num2str(n_nodes),'x',num2str(n_nodes),' Grid']])
ylabel('y [m]')
xlabel('x [m]')
colorbar

%% Explicit Euler routine
function T = explicit_euler(n_nodes, domain, boundary_conditions, solver)
    %unpack variables for readability
    %domain variables
    L = domain.L;
    %boundary_conditions variables
    T_x0y = boundary_conditions.T_x0y; %left
    dTdx_xLy = boundary_conditions.dTdx_xLy; %right
    T_xy0 = boundary_conditions.T_xy0; %bottom
    T_xyL = boundary_conditions.T_xyL; %top
    T_0internal = boundary_conditions.T_0internal; %internal points
    %solution control variables
    dt = solver.dt;
    t_stop = solver.t_stop;
    n_iterations = solver.n_iterations; %figure out how many time iterations to do

    T = T_0internal.*ones(n_nodes); %uniform (square) grid 
    %apply BCs to t0 entry
    T(1,:) = T_x0y; %left BC
    T(:,end) = T_xyL; %top BC
    T(:,1) = T_xy0; %bottom BC
    dx = L/n_nodes;
    dy = L/n_nodes;
    for n = 1:n_iterations      %loop over all timesteps
        for i = 2:n_nodes       %loop over 2:n x
            for j = 2:n_nodes-1 %loop over 2:n-1 y
                %compute x derivitive term
                if i == n_nodes %check if on right boundary
                    x_der = (dTdx_xLy - (T(i,j)...
                             - T(i-1,j)))/(dx^2); %apply right BC
                else
                    x_der = (1/(dx^2))*(T(i+1,j) - 2*T(i,j) + T(i-1,j));
                end
                %compute y derivitive term
                y_der = (1/(dy^2))*(T(i,j+1) - 2*T(i,j) + T(i,j-1));
                %compute phi @ t_n+dt
                T(i,j) = T(i,j) + (x_der + y_der)*dt;
            end %j loop
        end %i loop
    end %time loop
end % function explicit_euler