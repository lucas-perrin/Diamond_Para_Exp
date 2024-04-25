clear all
close all

%% DIAMOND STRATEGY WAVE EQ

%% ----- Space domain discretization

dx    = 0.01;
xspan = 0:dx:1;
mx    = length(xspan);

%% ----- State trajectory : continuous in time solution

alpha = 1;

x     = @(x,t) sin(alpha*4*pi.*x).*cos(4*pi.*t).*(x<1);
xdt   = @(x,t) -sin(alpha*4*pi.*x).*alpha*4*pi.*sin(4*pi.*t).*(x<1);
xdtt  = @(x,t) -sin(alpha*4*pi.*x).*(alpha^2)*16*pi^2.*cos(4*pi.*t).*(x<1);

%% ----- Semi-continuous discretization

D = ((mx-1)^2).*tridiag(1,-2,1,mx); D(1:1) = 1; D(end,end) = 1;
D(1,2) = 0; D(2,1) = 0; D(end,end-1) = 0; D(end-1,end) = 0;

A = [zeros(mx,mx) eye(mx); D zeros(mx,mx)];

X    = @(t) x(xspan',t);
Xdt  = @(t) xdt(xspan',t);
Xdtt = @(t) xdtt(xspan',t);

U    = @(t) [X(t)   ; Xdt(t) ];
Udt  = @(t) [Xdt(t) ; Xdtt(t)];

Xdxx = @(t) D*X(t);

G = @(t) [zeros(mx,1) ; (alpha^2).*Xdtt(t)-Xdxx(t)];

u0 = U(0);

%% ----- Recurrence for finding best dt and nested grids
if 0
%-----

% i = 1

dt_1      = 1e-1;
Tf        = 10;

error_end = [];
error_tot = [];
dt_s      = [];

for i = 5:12
    dt_i = dt_1/(2^i);
    tspan_i = 0:dt_i:10;
    [~,U_RK4_i] = odeRK4_inhom_ufunc(A,G,tspan_i,u0);
    error_end   = [error_end, norm(U_RK4_i(:,end) - U(tspan_i(:,end)))];
    error_tot   = [error_tot, sum(vecnorm(U_RK4_i - U(tspan_i)))];
    dt_s        = [dt_s, dt_i];
end

deb = find(~isnan(error_end), 1);

figure(1)
title('Error convergence of RK4 method (mx = 21)')
semilogy(dt_s(deb:end),error_end(deb:end),'b',...
    dt_s(deb:end),error_tot(deb:end),'r',...
    dt_s(deb:end),dt_s(deb:end).^3,'--r',...
    dt_s(deb:end),dt_s(deb:end).^4,'--b')
legend('error end','error total','dt^3','dt^4')
annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', sprintf("On prend : dt_min = 1e-1/(2^6) and dt_max = 1e-1/(2^10)"))

grid on

%-----
end
%% ----- Nested grids

% We take : dt_min = 1e-1/(2^5) and dt_max = 1e-1/(2^11)

%% ----- Observer setting

m_obs = floor(0.5*mx);
C0    = [zeros(m_obs,floor((mx-m_obs)/2)) eye(m_obs) zeros(m_obs,floor((mx-m_obs+1)/2))];
C     = [zeros(m_obs,mx) C0];
gain  = 200;
L     = gain * C';

Y     = @(t) C * U(t);

B     = eye(2*mx);

M     = A - L * C;
Gobs  = @(t) B * G(t) + L * Y(t);

xh0   = (40*rand(1)-20).*(xspan'+(2*rand(1)-1)).*(xspan'+(2*rand(1)-1)).*(xspan'+(2*rand(1)-1)).*(xspan'-1).*xspan'.*(0<xspan').*(xspan'<1);
uh0   = [xh0 ; xh0];

[V,E] = eig(M);

eigenvalues_M = diag(E);

mu = min(abs(eigenvalues_M(eigenvalues_M < 0)));
gamma = cond(V);

theorical_error = @(t) gamma.*exp(-t.*mu);

%% ----- gain and mu / number of observation relation
if 0
%-----

gain_list    = 2.^(-10:1:40);
m_obs        = floor(0.5*mx);
C0           = [zeros(m_obs,floor((mx-m_obs)/2)) eye(m_obs) zeros(m_obs,floor((mx-m_obs+1)/2))];
C            = [zeros(m_obs,mx) C0];
mu_gain_spd  = [];

for gg = gain_list
    L             = gg * C';
    M             = A - L * C;
    [V,E]         = eig(M);
    eigenvalues_M = diag(E);
    mu            = min(abs(eigenvalues_M(eigenvalues_M < 0)));
    mu_gain_spd   = [mu_gain_spd mu];
end

C           = [C0 zeros(m_obs,mx)];
mu_gain_pos = [];

for gg = gain_list
    L             = gg * C';
    M             = A - L * C;
    [V,E]         = eig(M);
    eigenvalues_M = diag(E);
    mu            = min(abs(eigenvalues_M(eigenvalues_M < 0)));
    mu_gain_pos   = [mu_gain_pos mu];
end

figure(2)
loglog(gain_list,mu_gain_spd,'r',gain_list,mu_gain_pos,'--b')
legend('speed observed','position observed')
grid on
xlabel('gain')
ylabel('mu')
title('Observer for wave equation : gain vs mu, 50% observed')

size_obs = [0.01 0.02 0.05 0.1 0.25 0.5 0.75 0.9 1];
gain     = 500;
mu_size_spd  = [];

for ss_obs = size_obs
    m_obs         = floor(ss_obs*mx);
    C0            = [zeros(m_obs,floor((mx-m_obs)/2)) eye(m_obs) zeros(m_obs,floor((mx-m_obs+1)/2))];
    C             = [zeros(m_obs,mx) C0];
    L             = gain * C';
    M             = A - L * C;
    [V,E]         = eig(M);
    eigenvalues_M = diag(E);
    mu            = min(abs(eigenvalues_M(eigenvalues_M < 0)));
    mu_size_spd   = [mu_size_spd mu];
end

mu_size_pos  = [];

for ss_obs = size_obs
    m_obs         = floor(ss_obs*mx);
    C0            = [zeros(m_obs,floor((mx-m_obs)/2)) eye(m_obs) zeros(m_obs,floor((mx-m_obs+1)/2))];
    C             = [C0 zeros(m_obs,mx)];
    L             = gain * C';
    M             = A - L * C;
    [V,E]         = eig(M);
    eigenvalues_M = diag(E);
    mu            = min(abs(eigenvalues_M(eigenvalues_M < 0)));
    mu_size_pos   = [mu_size_pos mu];
end

figure(3)
semilogy(size_obs,mu_size_spd,'r',size_obs,mu_size_pos,'--b')
legend('speed observed','position observed')
grid on
xlabel('size obs')
ylabel('mu')
title('Observer for wave equation : size obs vs mu, gain = 500')

%-----
end
%% ----- Luenberger for PDE
if 0
%-----

Tf        = 20;
dt_min    = (1e-1)/(2^5);
tspan_min = 0:dt_min:Tf;

[~,U_hat_min] = odeRK4_inhom_ufunc(M,Gobs,tspan_min,uh0);

dt_max    = (1e-1)/(2^11);
tspan_max = 0:dt_max:Tf;

[~,U_hat_max] = odeRK4_inhom_ufunc(M,Gobs,tspan_max,uh0);

figure(4)
semilogy(tspan_min,vecnorm(U_hat_min - U(tspan_min)),...
    tspan_max,vecnorm(U_hat_max - U(tspan_max)),...
    tspan_min(1:floor(length(tspan_min)/2)),theorical_error(tspan_min(1:floor(length(tspan_min)/2))),'--k')
legend('dt = (1e-1)/(2^5)','dt = (1e-1)/(2^11)','theorical_error','Interpreter','LateX')
title('Wave equation, || U hat - U ||, gain = 200, taille obs = 50%')
grid on

%-----
end
%% ----- Brouillon
if 0
%-----

tspan_1 = 1:0.001:10;

U_mat_1 = U(tspan_1);

[~,U_rk4_1] = odeRK4_inhom_ufunc(A,G,tspan_1,u0);

figure(10)
plot(tspan_1,vecnorm(U_mat_1 - U_rk4_1))

tspan_2 = 1:0.0001:10;

U_mat_2 = U(tspan_2);

[~,U_rk4_2] = odeRK4_inhom_ufunc(A,G,tspan_2,u0);

figure(10)
plot(tspan_2,vecnorm(U_mat_2 - U_rk4_2))

sum(vecnorm(U_mat_1 - U_rk4_1))
sum(vecnorm(U_mat_2 - U_rk4_2))

%-----
end