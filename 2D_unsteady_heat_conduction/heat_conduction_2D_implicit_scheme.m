% Solution of 2-D Unsteady Heat Conduction equation using Implicit Euler
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

%% calling Implicit Euler, plotting results
figure(1)
figure(2)
for k = 1:numel(grid_sizes) %loop over all grids
    n_nodes = grid_sizes(k);
    %make grid coords. for plotting
    x = linspace(0,domain.L,n_nodes);
    y = linspace(0,domain.L,n_nodes);
    %call implicit euler solver for current n_nodes, etc.
    T = implicit_euler(n_nodes, domain, boundary_conditions, solver);
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

figure(1)
title(['T(x=0.5,y) @ t=',num2str(solver.t_stop),'s', ', dt=',num2str(solver.dt)])
legend('location','Best')
box on, grid on
ylabel('T')
xlabel('y coordinate')

figure(2)
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

%% Implicit Euler routine
function T = implicit_euler(n_nodes, domain, boundary_conditions, solver)
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
    %var init.
    T = T_0internal.*ones(n_nodes); %uniform (square) grid 
    %apply BCs to t0 entry
    T(1,:) = T_x0y; %left BC
    T(:,end) = T_xyL; %top BC
    T(:,1) = T_xy0; %bottom BC
    dx = L/n_nodes;
    dy = L/n_nodes;
    for n = 1:n_iterations
        %ADI step 1 - sweeping over Y
        T_star = zeros(size(T));
        T_star(:,1) = T(:,1); %keep bottom BC.
        for j_level = 2:n_nodes-1 %iterate over j-levels
            %step 1 of ADI
            %make A1 matrix
            A1 = eye(n_nodes);
            %populate inner rows of A matrix
            for i = 2:n_nodes %x1->xn 
                if i == n_nodes %check if on right BC
                    AW = -1/(dx^2);
                    AP = (2/dt) + (1/(dx^2));
                    A1(i,i-1) = AW;
                    A1(i,i) = AP;
                else
                    AW = -1/(dx^2);
                    AP = (2/dt) + (2/(dx^2));
                    AE = -1/(dx^2);
                    %fill matrix
                    A1(i,i-1) = AW;
                    A1(i,i) = AP;
                    A1(i,i+1) = AE;
                end
            end
            %make Q vector
            Q1 = zeros(size(A1,1),1);
            Q1(1) = T(1,j_level); %apply left BC
            for ii = 2:n_nodes %loop over inner (and last) x nodes
                y_der = (T(ii,j_level+1)-2*T(ii,j_level)+T(ii,j_level-1))/(dy^2);
                Q1(ii) = (2/dt)*T(ii,j_level) + y_der;
            end
            %use TDMA routine to solve for T_star vector
            T_star(:,j_level) = TDMA_solve(A1,Q1);
            T_star(:,end) = T(:,end); %keep top BC.
            clear Q1 %clean up
        end %j-level loop
            
        %ADI step 2 - sweeping over X
        T_np1 = zeros(size(T_star));
        T_np1(1,:) = T_star(1,:); %keep left BC.
        for i_level = 2:n_nodes %iterate over i-levels
            %make A2 matrix
            %A2 = make_ADI_mat_step2(dy,dt,n_nodes);
            A2 = eye(n_nodes);
            for i = 2:n_nodes-1 %inner rows of matrix
                AW = -1/(dy^2);
                AP = (2/dt) + (2/(dy^2));
                AE = -1/(dy^2);
                %fill matrix
                A2(i,i-1) = AW;
                A2(i,i) = AP;
                A2(i,i+1) = AE;
            end
            %make Q vector
            Q2 = zeros(size(A2,1),1);
            Q2(1) = T_xy0; %apply bottom BC
            Q2(end) = T_xyL; %apply top BC
            for jj = 2:n_nodes-1
                if i_level == n_nodes %check if this i_level is right face
                    x_der = (dTdx_xLy - (T_star(i_level,jj)...
                             - T_star(i_level-1,jj)))/(dx^2); %apply right dTdx der. BC
                else
                    x_der = (T_star(i_level+1,jj)...
                             - 2*T_star(i_level,jj)...
                             + T_star(i_level-1,jj)...
                             )/(dx^2);
                end
                Q2(jj) = (2/dt)*T_star(i_level,jj) + x_der;
            end
            
            %use TDMA routine to solve for T_n+1 vector
            T_np1(i_level,:) = TDMA_solve(A2,Q2);
            clear Q2
        end %i-level loop
        T = T_np1; %update 'T' for next iteration
        clear T_np1 T_star %clean up
    end %time loop
end % function implicit_euler

%% numerical methods functions
function phi = TDMA_solve(A,Q) %Tri-diagonal Solver
    phi = zeros(size(Q));
    n = length(Q);
    %forward elimination
    for i = 2:n
        A(i,i) = A(i,i) - (A(i,i-1)/A(i-1,i-1))*A(i-1,i);
        Q(i) = Q(i) - (A(i,i-1)/A(i-1,i-1))*Q(i-1);
    end
    %back substitution
    phi(end) = Q(end)/A(n,n);
    for i = n-1:-1:1
        phi(i) = (Q(i) - A(i,i+1)*phi(i+1))/A(i,i);
    end
end