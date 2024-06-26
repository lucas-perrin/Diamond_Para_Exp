function [t_final,sol_final] = odeRK_nosave_inhom(order,A,U,tspan,x0,toler)

    if (nargin < 6)
        null_toler = true;
        toler = 1e-20;
    else
        null_toler = false;
    end

    n = length(tspan);

    sol_final_t = x0;

    if order == 4

        for t = 1:n-1
            
            dt = abs(tspan(t+1) - tspan(t));
            
            dt2 = dt/2;
            
            xk1 = A * sol_final_t + U(tspan(t));
            xk2 = A * (sol_final_t + dt2 * xk1) + U(tspan(t) + dt2);
            xk3 = A * (sol_final_t + dt2 * xk2) + U(tspan(t) + dt2);
            xk4 = A * (sol_final_t +  dt * xk3) + U(tspan(t+1));
            
            sol_final_tp1 = sol_final_t + (dt/6).*(xk1 + 2*xk2 + 2*xk3 + xk4);
                
            if (~ null_toler)*(norm(sol_final_tp1 - sol_final_t) < toler)
                sol_final_t = sol_final_tp1;
                t_final = tspan(t+1);
                break
            end
            
            t_final = tspan(t+1);
            sol_final_t = sol_final_tp1;

        end

    elseif order == 2

        for t = 1:n-1
        
            dt = abs(tspan(t+1) - tspan(t));
            
            xk1 = A * sol_final_t + U(tspan(t));
            xk2 = A * (sol_final_t + dt * xk1) + U(tspan(t) + dt);
            
            sol_final_tp1 = sol_final_t + (dt/2).*(xk1 + xk2);
                
            if (~ null_toler)*(norm(sol_final_tp1 - sol_final_t) < toler)
                sol_final_t = sol_final_tp1;
                t_final = tspan(t+1);
                break
            end
            
            t_final = tspan(t+1);
            sol_final_t = sol_final_tp1;

        end

        sol_final = sol_final_t;
    end

end