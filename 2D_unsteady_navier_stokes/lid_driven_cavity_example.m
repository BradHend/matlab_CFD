%lid-driven cavity example that uses 'unsteady_NS_2D' solver

clear, clc, close all

%bool flag for printing of convergence assessments
% print_convergence = true;
print_convergence = false;

%grid sizes to evaluate
% grid_sizes = [21, 41, 61]; %takes a bit of time
grid_sizes = [21]; %takes less time

%% Boundary Conditions for lid-driven cavity problem
%BCs for u velocity
boundary_conditions.u_x0y = 0.0; %left
boundary_conditions.u_xLy = 0.0; %right
boundary_conditions.u_xy0 = 0.0; %bottom
boundary_conditions.u_xyL = 1.0; %top
boundary_conditions.u_0internal = 0.0; %internal points
%BCs for v velocity
boundary_conditions.v_x0y = 0.0; %left
boundary_conditions.v_xLy = 0.0; %right
boundary_conditions.v_xy0 = 0.0; %bottom
boundary_conditions.v_xyL = 0.0; %top
boundary_conditions.v_0internal = 0.0; %internal points
%% Set up input params. for solver routine
%domain params.
domain.L = 1.0;
domain.nu = 0.01;
domain.rho = 1.0;
domain.mu = domain.rho*domain.nu;
%solver control
% dt = min([0.25*dx*dx/nu, 4.0*nu/u_top/u_top]);
solver.dt = 1e-3;
solver.print_convergence = print_convergence;
solver.steady_state_tol = 1e-5;

%% Calling unsteady_NS_2D solver given BCs, grid, solver params
for i = 1:numel(grid_sizes) %loop over grid size list
    n_nodes = grid_sizes(i);
    domain.n_nodes = n_nodes; %overide grid size
    
    %call solver with current grid size
    [u, v, P, plot_data] = unsteady_NS_2D(domain, boundary_conditions, solver);
    
    %remove ghost cells from plotting data
    u_plotting = u(2:end-1,2:end-1);
    v_plotting = v(2:end-1,2:end-1);
    P_plotting = P(2:end-1,2:end-1);
    %make grid coords. for plotting
    x = linspace(0, domain.L, grid_sizes(i));
    y = linspace(0, domain.L, grid_sizes(i));
    
    plot_idx = round(n_nodes/2); %find idx that corresponds to center of domain
    figure(1)
    hold on
    plot(u_plotting(plot_idx,:), y, 'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])
    
    figure(2)
    hold on
    plot(x, u_plotting(:,plot_idx), 'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])

    figure(3)
    hold on
    plot(v_plotting(plot_idx,:), y, 'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])
    
    figure(4)
    hold on
    plot(x, v_plotting(:,plot_idx), 'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])
    
    figure(5)
    hold on
    plot(P_plotting(plot_idx,:), y, 'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])
    
    figure(6)
    hold on
    plot(x, P_plotting(:,plot_idx), 'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])
    
    figure(7) %plot u-vel resid
    hold on
    plot(plot_data.conv_iteration_number, plot_data.u_resididual, 'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])
    figure(8) %plot v-vel resid
    hold on
    plot(plot_data.conv_iteration_number, plot_data.u_resididual, 'DisplayName',[num2str(n_nodes),'x',num2str(n_nodes),' Grid'])
end %loop over grid size list

%% Plotting
%compute vorticity for plotting
dx = domain.L/domain.n_nodes;
dy = dx;
vorticity = zeros(size(u(2:end-1,2:end-1)));
for i = 2:size(u,1)-1
    for j = 2:size(u,2)-1
        vorticity(i-1,j-1) = (v(i,j) - v(i-1,j))/dx ...
                        - (u(i,j) - u(i-1,j))/dy;
    end
end

figure(1)
title(['u(x=0.5,y) at steady state', ', dt=',num2str(solver.dt)])
legend('location','Best')
box on, grid on
xlabel('u [m/s]')
ylabel('y coordinate')
xlim([-1,1])

figure(2)
title(['u(x,y=0.5) at steady state', ', dt=',num2str(solver.dt)])
legend('location','Best')
box on, grid on
ylabel('u [m/s]')
xlabel('x coordinate')
ylim([-1,1])

figure(3)
title(['v(x=0.5,y) at steady state', ', dt=',num2str(solver.dt)])
legend('location','Best')
box on, grid on
xlabel('v [m/s]')
ylabel('y coordinate')
xlim([-1,1])

figure(4)
title(['v(x,y=0.5) at steady state', ', dt=',num2str(solver.dt)])
legend('location','Best')
box on, grid on
ylabel('v [m/s]')
xlabel('x coordinate')
ylim([-1,1])

figure(5)
title(['Pressure(x=0.5,y) at steady state', ', dt=',num2str(solver.dt)])
legend('location','Best')
box on, grid on
xlabel('Pressure')
ylabel('y coordinate')

figure(6)
title(['Pressure(x,y=0.5) at steady state', ', dt=',num2str(solver.dt)])
legend('location','Best')
box on, grid on
ylabel('Pressure')
xlabel('x coordinate')

figure(7)
title('Residuals - U Velocity')
legend('location','Best')
box on, grid on
xlabel('Iteration #')
ylabel('Residual')

figure(8)
title('Residuals - V Velocity')
box on, grid on
legend('location','Best')
xlabel('Iteration #')
ylabel('Residual')

figure(9)
contourf(x, y, u_plotting')
title(['U Velocity Contour: ',[num2str(n_nodes),'x',num2str(n_nodes),' Grid']])
ylabel('y [m]')
xlabel('x [m]')
colorbar

figure(10)
contourf(x, y, v_plotting')
title(['V Velocity Contour: ',[num2str(n_nodes),'x',num2str(n_nodes),' Grid']])
ylabel('y [m]')
xlabel('x [m]')
colorbar

figure(11)
imagesc(x, y, P_plotting')
title(['Cell-Centered Pressure: ',[num2str(n_nodes),'x',num2str(n_nodes),' Grid']])
ylabel('y [m]')
xlabel('x [m]')
set(gca, 'YDir', 'normal')
colorbar

% time = solver.dt.*[1:numel(plot_data.p_corner_BL)];
time = plot_data.p_time;
figure(12)
plot(time, plot_data.p_corner_BL)
title(['Pressure vs time, Bottom Left Cell : ',[num2str(n_nodes),'x',num2str(n_nodes),' Grid']])
ylabel('Pressure')
xlabel('time [s]')

figure(13)
plot(time, plot_data.p_corner_TL)
title(['Pressure vs time, Top Left Cell: ',[num2str(n_nodes),'x',num2str(n_nodes),' Grid']])
ylabel('Pressure')
xlabel('time [s]')

figure(14)
plot(time, plot_data.p_corner_TR)
title(['Pressure vs time, Top Right Cell: ',[num2str(n_nodes),'x',num2str(n_nodes),' Grid']])
ylabel('Pressure')
xlabel('time [s]')

figure(15)
plot(time, plot_data.p_corner_BR)
title(['Pressure vs time, Bottom Right Cell: ',[num2str(n_nodes),'x',num2str(n_nodes),' Grid']])
ylabel('Pressure')
xlabel('time [s]')

figure(16)
contourf(x, y, vorticity')
title(['Vorticity Contour: ',[num2str(n_nodes),'x',num2str(n_nodes),' Grid']])
ylabel('y [m]')
xlabel('x [m]')
cbar = colorbar;
cbar.Label.String = 'Vorticity';