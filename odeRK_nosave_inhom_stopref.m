function [t_final,sol_final,time] = odeRK_nosave_inhom_stopref(order,A,U,tspan,x0,toler,ref)

    time=0;

    tic

    n = length(tspan);
    sol_final = x0;
    
    time=time+toc;
    if order == 4

        for t = 1:n-1
            
            tic

            dt = abs(tspan(t+1) - tspan(t));
            
            dt2 = dt/2;
            
            xk1 = A * sol_final + U(tspan(t));
            xk2 = A * (sol_final + dt2 * xk1) + U(tspan(t) + dt2);
            xk3 = A * (sol_final + dt2 * xk2) + U(tspan(t) + dt2);
            xk4 = A * (sol_final +  dt * xk3) + U(tspan(t+1));
            
            sol_final = sol_final + (dt/6).*(xk1 + 2*xk2 + 2*xk3 + xk4);
            
            time=time+toc;

            if (norm(sol_final - ref(tspan(t+1))) < toler)

                tic

                t_final = tspan(t+1);

                time=time+toc;

                break
            end
            
            tic

            t_final = tspan(t+1);

            time=time+toc;
        end
    
    elseif order == 2
    
        for t = 1:n-1
            
            tic

            dt = abs(tspan(t+1) - tspan(t));
            
            xk1 = A * sol_final + U(tspan(t));
            xk2 = A * (sol_final + dt * xk1) + U(tspan(t) + dt);
            
            sol_final = sol_final + (dt/2).*(xk1 + xk2);
            
            time=time+toc;

            if (norm(sol_final - ref(tspan(t+1))) < toler)
                
                tic

                t_final = tspan(t+1);

                time=time+toc;
                break
            end
            
            tic

            t_final = tspan(t+1);
            
            time=time+toc;
        end
    
    end
    
end