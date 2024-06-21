function [tspan_bk,sol] = odeRK4_hom(A,tspan,x0)

m = length(x0);
n = length(tspan);

tau = tspan(2) - tspan(1); 

P4 = eye(m) + (tau.*A) + (1/2).*(tau.*A)^2 + (1/6).*(tau.*A)^3 + (1/24).*(tau.*A)^4;

sol = zeros(m,n);
sol(:,1) = x0;


for i = 2:n
    
    sol(:,i) = P4 * sol(:,i-1);
    
end

tspan_bk = tspan;

end
