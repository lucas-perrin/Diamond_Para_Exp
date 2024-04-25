function [tspan_bk,sol] = odeRK4_inhom_ufunc(A,U,tspan,x0,toler)

if (nargin < 5)
    null_toler = true;
    toler = 1e-20;
else
    null_toler = false;
end

m = length(x0);
n = length(tspan);

sol = zeros(m,n);

sol(:,1) = x0;

tspan_bk = tspan;

for t = 1:n-1
    
    dt = abs(tspan(t+1) - tspan(t));
    
    dt2 = dt/2;
    
    xk1 = A * sol(:,t) + U(tspan(t));
    xk2 = A * (sol(:,t) + dt2 * xk1) + U(tspan(t) + dt2);
    xk3 = A * (sol(:,t) + dt2 * xk2) + U(tspan(t) + dt2);
    xk4 = A * (sol(:,t) +  dt * xk3) + U(tspan(t+1));
    
    sol(:,t+1) = sol(:,t) + (dt/6).*(xk1 + 2*xk2 + 2*xk3 + xk4);
        
    if (~ null_toler)*(norm(sol(:,t+1) - sol(:,t)) < toler)
        tspan_bk = tspan(1:t+1);
        break
    end
    
    tspan_bk = tspan;
    
end

end