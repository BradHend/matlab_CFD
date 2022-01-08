% Solution of 1-D convection-Diffusion equation:
%   d_dx(rho*u*phi) = d_dx(gamma * dphi_dx)

clear, clc, close all

%% User configurable params.
n_nodes = 11;
L = 1.0;
rho = 1.0;
u = 1.0;

%BCs
phi0 = 0;
phiL = 1;
Pe = 50.0; %Pe to use for numerical solution
Pe_analytical = 50.0; %Pe to use for 'exact' solution
Gamma = rho*u*L/Pe;

%create grid
x0 = 0.0;
xL = L;
x = uniform_grid(x0, xL, n_nodes);

%% analytical solution for comparison
c = (phiL - phi0)/(exp(Pe_analytical)-1);
phi_analytical = @(x) phi0 + c*(exp(Pe_analytical*x./L)-1);

%% Example with Upwind differencing for convective term
%make A_c matrix, using UDS
A_c_UDS = make_A_c_UDS_matrix(x,rho,u);
%make A_d matrix, using CDS
A_d_CDS = make_A_d_CDS_matrix(x,Gamma);
%add convective and diffusive terms
A_full = A_c_UDS + A_d_CDS;
%fix phi0 and phiL coefficients
A_full(1,1) = 1;
A_full(end,end) = 1;

%make Q vector
Q = zeros(size(A_full,1),1);
Q(1) = phi0;
Q(end) = phiL;

%use TDMA routine to solve for [q1,q2,..., qn]
phi_solved_UDS = TDMA_solve(A_full,Q);

%% Example with CDS for convection term
%make A_c matrix, using CDS
A_c_CDS = make_A_c_CDS_matrix(x,rho,u);
%make A_d matrix, using CDS
A_d_CDS = make_A_d_CDS_matrix(x,Gamma);
%add convective and diffusive terms
A_full = A_c_CDS + A_d_CDS;
%fix phi0 and phiL coefficients
A_full(1,1) = 1;
A_full(end,end) = 1;

%make Q vector
Q = zeros(size(A_full,1),1);
Q(1) = phi0;
Q(end) = phiL;

%use TDMA routine to solve for [q1,q2,..., qn]
phi_solved_CDS = TDMA_solve(A_full,Q);

%plot results
figure
fplot(phi_analytical, [0,1], 'k')
hold on
plot(x, phi_solved_UDS, 'b')
hold on
plot(x, phi_solved_CDS, 'r')
legend({['Exact Solution: Pe=', num2str(Pe_analytical)], ...
    ['Finite Difference: Pe=', num2str(Pe),' (UDS Convective)'], ...
    ['Finite Difference: Pe=', num2str(Pe),' (CDS Convective)']},'location','Best')
title(['1D Convection-Diffusion: ',num2str(n_nodes), ' Point Grid'])
xlabel('x [m]')
ylabel('\phi')


%% Grid resolution study
clear phi_analytical
grid_sizes = 8:2:2^6;
Pe_grid_study = [18.0, 50.0];
for i = 1:numel(Pe_grid_study)
    %need to recalc. based on local Pe
    Pe_local = Pe_grid_study(i);
    Gamma = rho*u*L/Pe_local;
    c = (phiL - phi0)/(exp(Pe_local)-1);
    phi_analytical = @(x) phi0 + c*(exp(Pe_local*x./L)-1);
    for j = 1:numel(grid_sizes) %loop over all grid sizes in the study
        x = uniform_grid(x0, xL, grid_sizes(j));
        A_c_CDS = make_A_c_CDS_matrix(x,rho,u);
        %make A_d matrix, assuming CDS
        A_d_CDS = make_A_d_CDS_matrix(x,Gamma);
        %add convective and diffusive terms
        A_full = A_c_CDS + A_d_CDS;
        %fix phi0 and phiL coefficients
        A_full(1,1) = 1;
        A_full(end,end) = 1;
        %make Q vector
        Q = zeros(size(A_full,1),1);
        Q(1) = phi0;
        Q(end) = phiL;
        %use Gauss Elimination routine to solve for [q1,q2,..., qn]
        phi_solved_grid = TDMA_solve(A_full,Q);
        %compute error term
        err(i,j) = sum(abs(phi_analytical(x) - phi_solved_grid'));
    end
end
figure
plot(grid_sizes, err(1,:), grid_sizes, err(2,:))
title('Solution Error vs. Grid Points')
legend({['Pe=',num2str(Pe_grid_study(1)),' (CDS Convection)'],['Pe=',num2str(Pe_grid_study(2)),' (CDS Convection)']})
xlabel('Number of Grid Points')
ylabel('Error = |Exact - FiniteDiff.|')

%% generalized function for creating grid
function x = uniform_grid(x0, xL, n)
    x = linspace(x0,xL,n);
end

%% functions to form "A" matrix terms
function A_c = make_A_c_UDS_matrix(x,rho,u) %convection term, Upwind differencing scheme
    A_c = eye(length(x));
    for i = 2:length(x)-1 %inner rows of matrix
        AW = -max(rho*u,0)/(x(i)-x(i-1));
        AE = min(rho*u,0)/(x(i+1)-x(i));
        AP = -(AE+AW);
        %fill matrix
        A_c(i,i-1) = AW;
        A_c(i,i) = AP;
        A_c(i,i+1) = AE;
    end
end

function A_c = make_A_c_CDS_matrix(x,rho,u) %convection term, CDS
    A_c = eye(length(x));
    for i = 2:length(x)-1 %inner rows of matrix
        AW = -rho*u/(x(i+1)-x(i-1));
        AE = rho*u/(x(i+1)-x(i-1));
        AP = -(AE+AW);
        %fill matrix
        A_c(i,i-1) = AW;
        A_c(i,i) = AP;
        A_c(i,i+1) = AE;
    end
end

function A_d = make_A_d_CDS_matrix(x,Gamma) %diffusion term, CDS
    A_d = eye(length(x));
    for i = 2:length(x)-1 %inner rows of matrix
        AW = -2*Gamma/((x(i+1)-x(i-1))*(x(i)-x(i-1)));
        AE = -2*Gamma/((x(i+1)-x(i-1))*(x(i+1)-x(i)));
        AP = -(AE + AW);
        %fill matrix
        A_d(i,i-1) = AW;
        A_d(i,i) = AP;
        A_d(i,i+1) = AE;
    end
end

%% supporting numerical methods functions
function phi = TDMA_solve(A,Q) %Tri-diagonal matrix solver
    %init vars.
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