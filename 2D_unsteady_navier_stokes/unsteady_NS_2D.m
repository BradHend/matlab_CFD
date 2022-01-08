function [u, v, P, plot_data] = unsteady_NS_2D(domain, boundary_conditions, solver)
%solver for Incompressible Navier-Stokes equations on 2D uniform grid
%   explicit Euler scheme for time
%   variables are stored on co-located grids

    %extract variables from inputs
    %domain values
    L = domain.L;
    n_nodes = domain.n_nodes;
    nu = domain.nu;
    rho = domain.rho;
    mu = domain.mu;
    x_max = n_nodes;
    y_max = n_nodes;
    dx = L/n_nodes;
    dy = dx;
    dz = 1.0; %assume unit-depth domain
    dx_inv = 1.0/dx;
    dy_inv = 1.0/dy;
    %BCs for u velocity
    u_x0y = boundary_conditions.u_x0y; %left
    u_xLy = boundary_conditions.u_xLy; %right
    u_xy0 = boundary_conditions.u_xy0; %bottom
    u_xyL = boundary_conditions.u_xyL; %top
    u_0internal = boundary_conditions.u_0internal; %internal points
    %BCs for v velocity
    v_x0y = boundary_conditions.v_x0y; %left
    v_xLy = boundary_conditions.v_xLy; %right
    v_xy0 = boundary_conditions.v_xy0; %bottom
    v_xyL = boundary_conditions.v_xyL; %top
    v_0internal = boundary_conditions.v_0internal; %internal points
    %solution control variables
    dt = solver.dt;
    print_convergence = solver.print_convergence;
    if isfield(solver,'steady_state_tol') %check if user specified
        steady_state_tol = solver.steady_state_tol;
    else
        steady_state_tol = 1e-5;
    end
    alpha = 0.5;
    
    %initialize variable arrays
    P = zeros(n_nodes+2);
    P_np1 = P;
    u = u_0internal.*ones(n_nodes+2); %initialize all grid points
    v = v_0internal.*ones(n_nodes+2); %initialize all grid points
    u_last = u;
    v_last = v;
    
    converged = false;
    iter_count = 0;
    iter_count_idx = 0; %track how many residual info have been saved
    p_output_count_idx = 0; %track how many corner pressures have been saved
    conv_iteration_number = [];
    while (~converged)
            iter_count = iter_count + 1;
            %compute ghost point values for u
            u(1,:) = 2*u_x0y - u(2,:); %left-1
            u(end,:) = 2*u_xLy - u(end-1,:); %right+1
            u(:,end) = 2*u_xyL - u(:,end-1); %top+1
            u(:,1) = 2*u_xy0 - u(:,2); %bottom-1
            %compute ghost point values for v
            v(1,:) = 2*v_x0y - v(2,:); %left-1
            v(end,:) = 2*v_xLy - v(end-1,:); %right+1
            v(:,end) = 2*v_xyL - v(:,end-1); %top+1
            v(:,1) = 2*v_xy0 - v(:,2); %bottom-1

            %ensure BCs get applied to u* and v* before starting
            u_star = u;
            v_star = v;
            for i = 2:x_max+1      %x-loop, non-ghost interior cells
                for j = 2:y_max+1  %y-loop, non-ghost interior cells
                    %compute face-centered u-velocities
                    u_w = (u(i-1,j) + u(i,j))/2;
                    u_e = (u(i+1,j) + u(i,j))/2;
                    u_n = (u(i,j+1) + u(i,j))/2;
                    u_s = (u(i,j-1) + u(i,j))/2;
                    %compute face-centered v-velocities
                    v_n = (v(i,j+1) + v(i,j))/2;
                    v_s = (v(i,j-1) + v(i,j))/2;
                    %x convection terms
                    conv_x_w = u_w^2; %w face
                    conv_x_e = u_e^2; %e face
                    conv_x_n = u_n*v_n; %n face
                    conv_x_s = u_s*v_s; %s face
                    conv_x = rho*((conv_x_e*dy - conv_x_w*dy) + (conv_x_n*dx - conv_x_s*dx));
                    %x diffusion terms
                    %compute face-centered derivitive terms
                    du_dx_e = (u(i+1,j) - u(i,j))/dx;
                    du_dx_w = (u(i,j) - u(i-1,j))/dx;
                    du_dy_n = (u(i,j+1) - u(i,j))/dy;
                    du_dy_s = (u(i,j) - u(i,j-1))/dy;
                    diff_x_w = dy*du_dx_w; %w face
                    diff_x_e = dy*du_dx_e; %e face
                    diff_x_n = dx*du_dy_n; %n face
                    diff_x_s = dx*du_dy_s; %s face
                    diff_x = (diff_x_e - diff_x_w) + (diff_x_n - diff_x_s);
                    %compute u*
                    u_star(i,j) = u(i,j) + (dt/(rho*dx*dy))*(-conv_x + mu*diff_x);
                end
            end
            for i = 2:x_max+1      %x-loop, non-ghost interior cells
                for j = 2:y_max+1  %y-loop, non-ghost interior cells
                    %compute face-centered u-velocities
                    u_w = (u(i-1,j) + u(i,j))/2;
                    u_e = (u(i+1,j) + u(i,j))/2;
                    %compute face-centered v-velocities
                    v_w = (v(i-1,j) + v(i,j))/2;
                    v_e = (v(i+1,j) + v(i,j))/2;
                    v_n = (v(i,j+1) + v(i,j))/2;
                    v_s = (v(i,j-1) + v(i,j))/2;

                    %y momentum terms
                    % compute face-centered derivitive terms
                    dv_dx_e = (v(i+1,j) - v(i,j))/dx;
                    dv_dx_w = (v(i,j) - v(i-1,j))/dx;
                    dv_dy_n = (v(i,j+1) - v(i,j))/dy;
                    dv_dy_s = (v(i,j) - v(i,j-1))/dy;
                    conv_y_w = u_w*v_w; %w face
                    conv_y_e = u_e*v_e; %e face
                    conv_y_n = v_n^2; %n face
                    conv_y_s = v_s^2; %s face
                    conv_y = rho*((conv_y_e*dy - conv_y_w*dy) + (conv_y_n*dx - conv_y_s*dx));
                    %y diffusion terms
                    diff_y_w = dy*dv_dx_w; %w face
                    diff_y_e = dy*dv_dx_e; %e face
                    diff_y_n = dx*dv_dy_n; %n face
                    diff_y_s = dx*dv_dy_s; %s face
                    diff_y = (diff_y_e - diff_y_w) + (diff_y_n - diff_y_s);
                    %compute v*
                    v_star(i,j) = v(i,j) + (dt/(rho*dx*dy))*(-conv_y + mu*diff_y);
                end
            end

            %solve for pressure term at next time step
            for i = 2:x_max+1 %loop over x-coords
                for j = 2:y_max+1 %loop over y-coords
                    term1 = (u_star(i+1,j) - u_star(i-1,j))/(2*dx);
                    term2 = (v_star(i,j+1) - v_star(i,j-1))/(2*dy);
                    Q(i-1,j-1) = (rho/dt)*(term1 + term2);
                end
            end
            %call SOR solver
            P_solved = solve_presure_poisson(Q, x_max, y_max, dx, dy);
            P_np1(2:end-1,2:end-1) = P_solved;
            %compute ghost point values for P that gives zero gradient
            P_np1(:,end) = P_np1(:,end-1); %top+1
            P_np1(:,1)   = P_np1(:,2); %bottom-1
            P_np1(1,:)   = P_np1(2,:); %left-1
            P_np1(end,:) = P_np1(end-1,:); %right+1

            for i = 2:x_max+1 %loop over x-coords "last"
                for j = 2:y_max+1 %loop over y-coords "first"
                    %compute correction for x-vel
                    dP_dx_ij = (P_np1(i+1,j) - P_np1(i-1,j))/(2*dx);
                    %compute correction for y-vel
                    dP_dy_ij = ((P_np1(i,j+1) - P_np1(i,j-1))/(2*dy));
                    %compute corrected velocities
                    u(i,j) = (1-alpha)*u(i,j) + alpha*(u_star(i,j) - (dt/rho)*dP_dx_ij);
                    v(i,j) = (1-alpha)*v(i,j) + alpha*(v_star(i,j) - (dt/rho)*dP_dy_ij);
                end
            end
            %check for convergence of u/v
            u_res = max(max(abs(u-u_last)));
            v_res = max(max(abs(v-v_last)));
            if (u_res < steady_state_tol) && (v_res < steady_state_tol)
                converged = true;
                break
            else
                %check if convergence data needs to be saved/printed
                if mod(iter_count, 100) == 0 %output residual info every 100 iterations
                    iter_count_idx = iter_count_idx + 1;
                    if print_convergence
                        disp(['u_redidual @iter ',num2str(iter_count), ' ', num2str(u_res)])
                        disp(['v_redidual @iter ',num2str(iter_count), ' ', num2str(v_res)])
                    end
                    u_resididual(iter_count_idx) = u_res;
                    v_resididual(iter_count_idx) = v_res;
                    conv_iteration_number(iter_count_idx) = iter_count;
                end
            end
            if mod(iter_count, 10) == 0 %output pressure info every 10 iterations
                p_output_count_idx = p_output_count_idx + 1;
                p_time(p_output_count_idx) = iter_count.*solver.dt; %save time info
                p_corner_BL(p_output_count_idx) = P(2,2);
                p_corner_TL(p_output_count_idx) = P(end-1,2);
                p_corner_BR(p_output_count_idx) = P(end-1,2);
                p_corner_TR(p_output_count_idx) = P(end-1,end-1);
            end
            
            %update solved flow variables for next timestep
            P = P_np1; 
            u_last = u;
            v_last = v;
    end %time loop
    
    %collect plot data for return
    plot_data.u_resididual = u_resididual;
    plot_data.v_resididual = v_resididual;
    plot_data.conv_iteration_number = conv_iteration_number;
    plot_data.p_time = p_time;
    plot_data.p_corner_BL = p_corner_BL;
    plot_data.p_corner_TL = p_corner_TL;
    plot_data.p_corner_BR = p_corner_BR;
    plot_data.p_corner_TR = p_corner_TR;
end %end unsteady_NS_2D

%% Poisson pressure solver
function p = solve_presure_poisson(Q,nx,ny,dx,dy)
    n_ghost = 1; %number of ghost cells to use
    p = zeros([ny,nx]);
    Ap = zeros([ny,nx]);
    Ae = (1.0/(dx*dx)).*ones([ny,nx]);
    As = (1.0/(dy*dy)).*ones([ny,nx]);
    An = (1.0/(dy*dy)).*ones([ny,nx]);
    Aw = (1.0/(dx*dx)).*ones([ny,nx]);
    %set left wall coefs
    Aw(1,:) = 0.0;
    %set right wall coefs
    Ae(end,:) = 0.0;
    %set top wall coefs
    An(:,end) = 0.0;
    %set bottom wall coefs
    As(:,1) = 0.0;
    Ap = -(Aw + Ae + An + As);
    alpha = 0.25;
    iter = 0;
    error = 1e2;
    tol = 1e-5;
    maxiter = 2000;
    while (error > tol) && (iter < maxiter) %loop until converged, or max iter
        pn = p; %make local copy for convergence check
        for i = 1:nx     %loop over x
            for j = 1:ny %loop over y
                %get coefficients for current cell
                ap = Ap(j,i);
                an = An(i,j);
                as = As(i,j);
                ae = Ae(i,j);
                aw = Aw(i,j);
                sum_term = 0.0; %var. init
                %check if on x-bounds
                if i == 1 %only add east face
                    sum_term = sum_term + ae*p(i+1,j);
                elseif i == nx %only add west face
                    sum_term = sum_term + aw*p(i-1,j);
                else %add both east and west
                    sum_term = sum_term + ae*p(i+1,j) + aw*p(i-1,j);
                end
                %check if on y-bounds
                if j == 1 %only add north face
                    sum_term = sum_term + an*p(i,j+1);
                elseif j == ny %only add south face
                    sum_term = sum_term + as*p(i,j-1);
                else %add both north and south
                    sum_term = sum_term + an*p(i,j+1) + as*p(i,j-1);
                end
                rhs = Q(i,j) + ((1-alpha)/alpha)*ap*p(i,j) - sum_term;
                %solve for next iteration of pressure term
                p(i,j) = (alpha/ap)*rhs;
            end
        end
        error = max(max(abs(p - pn)));
        iter = iter + 1;
    end
end