%===============================================================
% 
%
%
% 
% 
% 
% 
%===============================================================

clear all
close all

fontSize = 16;

%===============================================================
% if true, runs a RK4 scheme to see if the cheme matches the
% true solution and if the 
test_state_obs = true;

%===============================================================
% Observer parameters
fprintf('Setting parameters... \n')

x0       = 0;           % start of space domain
Nx       = 128;         % Number of points in space domain
L        = 1;           % length of space domain
dt       = 1e-3;        % timestep
size_obs = 1;           % value between 0 and 1
g        = 9.81;        % gravity constant
T_start  = 0;          % start time
T        = 10;           % time window
T_end    = T_start + T; % end time
N_s      = 10;          % Number of wave frequency in the solution wave

p        = 16;  % #CPU
gain     = 1;  % gamma parameter

%===============================================================
% Space Discretization
fprintf('Setting space grid... \n')

xspan = linspace(x0, x0+L, Nx); % space grid
dx    = xspan(2) -xspan(1);      % step in space
dxsp  = L/(Nx-1);
xsp   = 0:dxsp:L;

%===============================================================
% Time Discretization
fprintf('Setting time grid... \n')

tspan = T_start:dt:T_end;
Nt    = length(tspan);

%===============================================================
% Evolution Matrix
fprintf('Setting evolution matrix... \n')

A_dz = Matrix_A_dz(Nx,dxsp);
A  = [zeros(Nx), -g.*eye(Nx); A_dz, zeros(Nx)];

%===============================================================
% True Solution
fprintf('Setting true solution... \n')

a_s     = ones(1,N_s).*[N_s:-1:1]./(N_s^2)./4;
kw_s    = 2*pi*[1:N_s]/L;
omega_s = sqrt(kw_s.*g);

eta_s   = @(t) sum(a_s'.*sin(omega_s'.*t - kw_s'*xsp));
phi0_s  = @(t) sum(a_s'.*(g.*diag(1./omega_s)*cos(omega_s'.*t - kw_s'.*xsp)));

U_s     = @(t) [reshape(phi0_s(t),Nx,1);reshape(eta_s(t),Nx,1)];

G_null  = @(t) zeros(2*Nx,1);

eta_s_dt  = @(t) sum(a_s'.*(omega_s'.*cos(omega_s'.*t - kw_s'*xsp)));
phi0_s_dt = @(t) sum(-g.*a_s'.*(sin(omega_s'.*t - kw_s'*xsp)));
U_dt      = @(t) [reshape(phi0_s_dt(t),Nx,1);reshape(eta_s_dt(t),Nx,1)];
G         = @(t) U_dt(t) - A*U_s(t);

Ugrid_s = zeros(2*Nx,Nt);

for t=1:Nt
    Ugrid_s(:,t) = U_s(tspan(t));
end

U0 = U_s(T_start);

%===============================================================
% Observer Setting
fprintf('Setting observer... \n')

U0obs      = zeros(2*Nx,1);
%obs_map    = ones(1,Nx);
%[m_obs,Cc] = GetC(obs_map);
%C          = [zeros(m_obs,Nx) Cc];
C          = [zeros(Nx),eye(Nx)];
Lmat       = gain * C';
M          = A - Lmat*C;

Y_s        = @(t) C * U_s(t);
Gobs_s     = @(t) G(t) + Lmat*Y_s(t);

%===============================================================
% Tests
fprintf('Doing tests... \n')

if test_state_obs
    fprintf('  Doing test state... \n')
    fprintf('    Running scheme... \n')

    [~,Ustate_test]  = odeRK4_inhom_ufunc(sparse(A),G,tspan,U0);
    error_state_test = vecnorm(Ustate_test - Ugrid_s);

    fprintf('  Doing test obs... \n')
    fprintf('    Running scheme... \n')

    [~,Uobs_test] = odeRK4_inhom_ufunc(sparse(M),Gobs_s,tspan,zeros(2*Nx,1));
    error_obs_test   = vecnorm(Uobs_test - Ugrid_s);
    error_vec     = abs(Uobs_test(1:Nx,:) - Ugrid_s(1:Nx,:));

    fprintf('  Doing test error equation... \n')
    fprintf('    Running scheme... \n')

    [~,Uerr_test] = odeRK4_inhom_ufunc(sparse(M),G_null,tspan,U0);
    error_eq_test    = vecnorm(Uerr_test);

    fprintf('  Finding theorical convergence rate... \n')

    [V_m,E_m]      = eig(full(M));
    cond_m         = cond(V_m);
    zz             = abs(diag(E_m))>1e-10;
    E_m            = diag(E_m);
    E_m            = E_m(zz);
    mu             = min(abs(E_m));
    error_th_test  = @(t) cond_m*exp(-(t-T_start).*mu);

    fprintf('  Plotting... \n')
    
    figure(1)
    semilogy(tspan,error_obs_test,'r',...
    tspan,error_state_test,'b',...
    tspan,error_eq_test,'g',...
    tspan,error_th_test(tspan),'--k')
    legend('error observer scheme','error state scheme','error equation scheme','theorical error', 'Interpreter', 'latex')
    ylim([1e-12, 1e1])
    xlabel("time $t$", 'Interpreter', 'latex','FontSize', fontSize)
    ylabel("$|| \epsilon(t) ||$", 'Interpreter', 'latex','FontSize', fontSize,'rotation',0)
    title(['(LWWE) norm error convergence plot, $N_x = 128$, $\Delta_t = 10^{-5}$, $T_0 = 15$'], 'Interpreter', 'latex','FontSize', fontSize)

    figure(2)
    plot(xspan,Uerr_test(1:Nx,end))
end