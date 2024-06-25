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
% if true, runs a loop to find the upper and lower bounds of the
% nested grids 
test_nested_grids = true;
% if true, runs a test to see if the defined observer converges
% at the expected convergence rate
test_obs_conv = true;
% if true, runs a test to see if the basic ParaExp algorithm applied
% to the observer beahves as expected
test_ParaExp = false;
% if true, runs a test to see if the ParaExp alogorithm that only
% returns the end value of the interval applied to the observer 
% behaves as expected
test_ParaExp_end = false;
% if true, display a review of the current during the strategy
verbose = true;

%===============================================================
% Tunable parameters
fprintf('\nSetting parameters... \n')

% Gobal parameters

% order of the RK method, can be 2 or 4
order = 4;
% grid discretization, can be 0.1 or 0.1/2
h = 0.1/2;
% diffusive coefficient of the heat equation
alpha_heat = 0.05;
% radius of the circle inside the [0,1]^2 square that represents the
% observation domain
size_obs  = 0.25;
% observer gain
gain = 0;%2e2;

% Tests parameters

% number of computers for tests
N_tests = 4;
% final computation time for tests
Tf = 20;

% Strategy parameters

% final computation time strategy
Tf_strategy    = 20;
% number of computers
N_strategy_list = 2.^[1:8];
% stopping tolerance criteria
toler_strategy = 1e-13;

T_strategy_list = [0.2,0.5,1,2,5];

%===============================================================
% Space Discretization
fprintf('\nSetting space grid... \n')

xspanS = 0:h:1;
yspanS = 0:h:1;
xspan  = xspanS(2:end-1);
yspan  = yspanS(2:end-1);
mx     = length(xspan);
my     = length(yspan);

%===============================================================
% Evolution Matrix
fprintf('\nSetting evolution matrix... \n')

Id    = eye(my);
D2    = (mx*my)*tridiag(1,-4,1,my);
D1    = (mx*my)*eye(my);

D_diag    = D2;
D_diag_up = D1;

for i = 1:mx-2
    D_diag    = blkdiag(D_diag,D2);
    D_diag_up = blkdiag(D_diag_up,D1);
end
D_diag = blkdiag(D_diag,D2);

D_diag_down = flip(blkdiag(flip(D_diag_up),zeros(my,my))); %-1
D_diag_up   = flip(blkdiag(zeros(my,my),flip(D_diag_up))); %+1

A = alpha_heat * sparse(D_diag + D_diag_up + D_diag_down);

%===============================================================
% True Solution
fprintf('\nSetting true solution... \n')

x    = @(x,y,t)      sin(pi.*y')*sin(pi.*x).*exp(-10.*t) + (1/2).*atan(t)    .*(exp(-20.*(y'-1/2).^2)*exp(-20.*(x-1/2).^2)).*(sin(pi.*y')*sin(pi.*x));
xdt  = @(x,y,t) -10.*sin(pi.*y')*sin(pi.*x).*exp(-10.*t) + (1/(1+t^2)).*(1/2).*(exp(-20.*(y'-1/2).^2)*exp(-20.*(x-1/2).^2)).*(sin(pi.*y')*sin(pi.*x));

X    = @(t) x(xspan,yspan,t);
Xdt  = @(t) xdt(xspan,yspan,t);

Xdxx = @(t) reshape(A*reshape(X(t),mx*my,1),my,mx);
G    = @(t) Xdt(t)-Xdxx(t);
G    = @(t) reshape(G(t),mx*my,1);

X_reshape = @(t) reshape(X(t),mx*my,1);

x0   = X(0);

%===============================================================
% Tests
fprintf('\nDoing tests... \n')

if test_nested_grids
    fprintf('  Doing test nested grids... \n')

    dt_1      = 1e-1;

    error_end = [];
    error_tot = [];
    dt_s      = [];

    fprintf('    Doing loop... \n')

    for i = 0:11
        dt_i = dt_1/(2^i);
        fprintf("      i = %d \n",i)
        tspan_i = 0:dt_i:Tf;
        X_sol_i = Solvetspan(X,tspan_i);
        [~,X_RK4_i] = odeRK_inhom_ufunc(order,A,G,tspan_i,x0(:));
        error_end   = [error_end, norm(X_RK4_i(:,end)-X_sol_i(:,end))];
        error_tot   = [error_tot, sum(vecnorm(X_RK4_i - X_sol_i))];
        dt_s        = [dt_s, dt_i];
    end

    min_pow_grid = find(error_end < 1, 1);
    fprintf(['      min_pow_grid = ',num2str(min_pow_grid),'\n'])

    fprintf('    Plotting... \n')

    figure(1)
    loglog(dt_s(min_pow_grid:end),error_end(min_pow_grid:end),'-x',...
        dt_s(min_pow_grid:end),error_tot(min_pow_grid:end),'-x',...
        dt_s(min_pow_grid:end),dt_s(min_pow_grid:end).^2,'-',...
        dt_s(min_pow_grid:end),dt_s(min_pow_grid:end).^3,'-',...
        dt_s(min_pow_grid:end),dt_s(min_pow_grid:end).^4,'-')
    title('Error convergence of RK4 method (mx = 9)')
    legend('error end','error total','dt^2','dt^3','dt^4')
    annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', sprintf("On prend : dt_ min = 1e-1/(2^ 6) and dt_ max = 1e-1/(2^ 10)"))
    grid on
end

%===============================================================
% Nested grids setting
fprintf('\nSetting nested grids... \n')

if h==0.1
    if order==4
        min_pow_grid = 5;
        max_pow_grid = 10;
    elseif order==2
        min_pow_grid = 6;
        max_pow_grid = 16;
    end
elseif h==0.1/2
    if order==4
        min_pow_grid = 7;
        max_pow_grid = 10;
    elseif order==2
        min_pow_grid = 7;
        max_pow_grid = 16;
    end
end

%===============================================================
% Observer Setting
fprintf('\nSetting observer... \n')

%obs_map    = ((5/9) - (yspan - (1/2)).^2)'*((5/9) - (xspan - (1/2)).^2) >= (1/2)^2;    % big middle circle
%obs_map    = ((5/9) - (yspan - (1/2)).^2)'*((5/9) - (xspan - (1/2)).^2) >=(6/11)^2;    % small midle circle
%obs_map    = (7/11)*sin(2*pi.*yspan)'*(7/11)*sin(pi.*xspan) >= (6/11)^2;               % shifted

[XX, YY] = meshgrid(xspan,yspan);

obs_map = ((XX - 0.5).^2) + ((YY - 0.5).^2) < size_obs^2;

[m_obs,C]  = GetC(obs_map);
L          = gain * C';
 
Y     = @(t) C * reshape(X(t),mx*my,1);

B     = eye(mx*my);

M     = A - L * C;
Gobs  = @(t) B * reshape(G(t),mx*my,1) + L * Y(t);

[V,E] = eig(M);
mu    = min(abs(diag(E)));
gamma = cond(V);

x0obs = -x0;

theorical_error = @(t) gamma.*exp(-t.*mu).*norm(x0 - x0obs);

M = sparse(M);

fprintf('  Plotting... \n')
figure(2)
spy(obs_map)
title('Observer Map')

%===============================================================
% Tests
fprintf('\nDoing tests... \n')

if test_obs_conv
    fprintf("  Doing test observer convergence... \n")

    dt_min    = (1e-1)/(2^min_pow_grid);
    tspan_min = 0:dt_min:Tf;

    fprintf("    Computing serial with dt min = 1/(2^%i)... \n",min_pow_grid)

    [~,X_hat_min] = odeRK_inhom_ufunc(order,M,Gobs,tspan_min,x0obs(:));

    dt_max    = (1e-1)/(2^max_pow_grid);
    tspan_max = 0:dt_max:Tf;

    fprintf("    Computing serial with dt min = 1/(2^%i)... \n",max_pow_grid)

    [~,X_hat_max] = odeRK_inhom_ufunc(order,M,Gobs,tspan_max,x0obs(:));

    fprintf("    Plotting... \n")

    figure(3)
    semilogy(tspan_min,vecnorm(X_hat_min - Solvetspan(X,tspan_min)),...
        tspan_max,vecnorm(X_hat_max - Solvetspan(X,tspan_max)),...
        tspan_min(1:floor(length(tspan_min)/2)),theorical_error(tspan_min(1:floor(length(tspan_min)/2))),'--k')
    legend(['dt = (1e-1)/(2^',num2str(min_pow_grid),')'],['dt = (1e-1)/(2^',num2str(max_pow_grid),')'],'theorical_error','Interpreter','LateX')
    title(['Heat equation, || X hat - X ||, gain = ',num2str(gain),', taille obs = 50%'])
    %annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', sprintf("PARFAIT !"))
    grid on

end

if test_ParaExp
    fprintf("  Doing test ParaExp with p = %i... \n",N_tests)

    dt_max    = (1e-1)/(2^10);
    tspan_max = 0:dt_max:Tf;
    
    fprintf("    Computing ParaExp with dt max = 1/(2 10)... \n")

    [X_hat_para_max,time_para_test] = ParaExp(order,N_tests,tspan_max,M,Gobs,x0obs(:));

    fprintf("    Computing serial with dt max = 1/(2 10)... \n")

    tic;
    [~,X_hat_max] = odeRK_inhom_ufunc(order,M,Gobs,tspan_max,x0obs(:));
    serial_time_test = toc;

    fprintf(['      speedup = ',num2str(serial_time_test/serial_time_test),'\n'])

    figure(4)
    semilogy(tspan_max,vecnorm(X_hat_max - Solvetspan(X,tspan_max)),...
        tspan_max,vecnorm(X_hat_para_max - Solvetspan(X,tspan_max)),...
        tspan_max(1:floor(length(tspan_max)/2)),theorical_error(tspan_max(1:floor(length(tspan_max)/2))),'--k')
    legend('serial',['parallel p = ',num2str(N_tests)],'theorical error','Interpreter','LateX')
    title('|| X hat para - X ||, dt = (1e-1)/(2^ 10), typz 2 : expm')
    %annotation('textbox', [0.2, 0.2, 0.1, 0.1], 'String', sprintf("Comme avant : petite erreur à l'attaque de l'exponentielle \n peut-être raffiner le type 2 ?"))
    grid on
end

if test_ParaExp_end
    fprintf("  Doing test ParaExp_end with p = %i... \n",N_tests)

    dt_max    = (1e-1)/(2^max_pow_grid);
    tspan_max = 0:dt_max:Tf;
    
    fprintf("    Computing ParaExp with dt max = 1/(2^%i)... \n",max_pow_grid)

    [X_hat_para_max_end,time_para_end_test] = ParaExp_end(order,N_tests,tspan_max,M,Gobs,x0obs(:));

    fprintf("    Computing serial with dt max = 1/(2^%i)... \n",max_pow_grid)

    tic;
    [~,X_hat_max] = odeRK_inhom_ufunc(order,M,Gobs,tspan_max,x0obs(:));
    serial_time_test = toc;

    fprintf(['      error ParaExp_end at end = ',num2str(norm(X_hat_para_max_end -  Solvetspan(X,tspan_max(end)))),'\n'])
    fprintf(['      error Serial at end = ',num2str(norm(X_hat_para_max_end -  Solvetspan(X,tspan_max(end)))),'\n'])
    fprintf(['      speedup = ',num2str(serial_time_test/time_para_end_test),'\n'])
end

%===============================================================
% Data Assimilation Strategy
fprintf('\nData assimilation strategy... \n')
if 0

fprintf("  Doing serial scheme...\n")

dt_max    = (1e-1)/(2^max_pow_grid);
tspan_max = 0:dt_max:Tf_strategy;
nt_max    = length(tspan_max);

[tspan_serial,X_hat_max,time_serial] = odeRK_inhom_ufunc_stop(order,M,Gobs,tspan_max,x0obs(:),toler_strategy,X_reshape);

fprintf("    t_end serial = %d \n",tspan_serial(end))

% (dt_l)^4 <= exp(-mu*(l-1)*T) * (dt_1)^4

dt_list   = (1e-1)./(2.^[min_pow_grid:1:max_pow_grid]);

for ti = 1:length(T_strategy_list)
    
    T              = T_strategy_list(ti);
    length_windows = floor(Tf_strategy/T);
    dim_int_max    = nt_max/length_windows;
            
    for Ni = 1:length(N_strategy_list)
        
        N = N_strategy_list(Ni);
        %----- Serial scheme
        
        fprintf("  ==== T = %d || N = %d ==== \n",T,N)
        
        %----- Choix facile
        
        i=1;
        dt      = dt_list(i);
        tspan   = 0:dt:T;
        dim_int = dim_int_max;
        nt      = nt_max;
        x0obsp  = x0obs(:);
        
        fprintf("    Doing parallel scheme...\n")
        
        time_para = 0;
    
        %----- Iteration on the windows
        for l = 1:length_windows
            
            fprintf("      Window number = %d / %d \n",l,length_windows)
    
            if l == 1
                dt = dt_list(1);
                W_l = 0:dt:T;
            else
                if l < length_windows
                    cond_nested_grid = dt_list < ((dt^order)*exp(-mu))^(1/order);
                    if sum(cond_nested_grid)==0
                        dt = dt_list(end);
                    else
                        i  = find(1==cond_nested_grid,1);
                        dt = dt_list(i);
                    end
                    W_l = (l-1)*T:dt:l*T;
                else
                    W_l = l*T:dt:Tf;
                end
            end

            fprintf('        i=%i\n',i)
    
            %----- ParaExp Solver
    
            [X_hat_para_end_l,time_l] = ParaExp_end(order,Ni,W_l,M,Gobs,x0obsp);
            
            time_para = time_para + time_l; 
            
            if norm(X_hat_para_end_l - X_reshape(W_l(end))) < toler_strategy
                % tspan_para = tspan(1:(l<length_windows)*l*dim_int + (l==length_windows)*nt);
                tspan_para_end = W_l(end);
                break
            else
                % tspan_para = tspan;
                tspan_para_end = tspan(end);
            end
            
            x0obsp = X_hat_para_end_l;
            
        end
        %----- Las window
    
        if verbose
            fprintf("speedup = %d \n\n",(time_serial/time_para))

            fprintf("time serial = %d \n",time_serial)
            fprintf("time parall = %d \n",time_para)
            fprintf("endtime serial = %d \n",tspan_serial(end))
            fprintf("endtime parall = %d \n",tspan_para_end)
            fprintf("error serial = %d \n",norm(X_hat_max(end) - X_reshape(tspan_serial(end))))
            fprintf("error parall = %d \n",norm(X_hat_para_end_l - X_reshape(tspan_para_end)))
            
            [time_serial time_para tspan_serial(end) tspan_para_end]
            
            Results(ti,Ni) = (time_serial/time_para);

        end
    end

end

end