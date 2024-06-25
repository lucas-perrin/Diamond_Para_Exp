function [tspan_bk,sol,time] = odeRK_inhom_ufunc_stop(order,A,U,tspan,x0,toler,ref)

    time=0;

    tic

    m = length(x0);
    n = length(tspan);
    
    sol = zeros(m,n);
    
    sol(:,1) = x0;
    
    tspan_bk = tspan;

    time=time+toc;
    if order == 4

        for t = 1:n-1
            
            tic

            dt = abs(tspan(t+1) - tspan(t));
            
            dt2 = dt/2;
            
            xk1 = A * sol(:,t) + U(tspan(t));
            xk2 = A * (sol(:,t) + dt2 * xk1) + U(tspan(t) + dt2);
            xk3 = A * (sol(:,t) + dt2 * xk2) + U(tspan(t) + dt2);
            xk4 = A * (sol(:,t) +  dt * xk3) + U(tspan(t+1));
            
            sol(:,t+1) = sol(:,t) + (dt/6).*(xk1 + 2*xk2 + 2*xk3 + xk4);
            
            time=time+toc;

            if (norm(sol(:,t+1) - ref(tspan(t+1))) < toler)

                tic

                tspan_bk = tspan(1:t+1);
                sol = sol(:,1:t+1);

                time=time+toc;

                break
            end
            
            tic

            tspan_bk = tspan;

            time=time+toc;
        end
    
    elseif order == 2
    
        for t = 1:n-1
            
            tic

            dt = abs(tspan(t+1) - tspan(t));
            
            xk1 = A * sol(:,t) + U(tspan(t));
            xk2 = A * (sol(:,t) + dt * xk1) + U(tspan(t) + dt);
            
            sol(:,t+1) = sol(:,t) + (dt/2).*(xk1 + xk2);
            
            time=time+toc;

            if (norm(sol(:,t+1) - ref(tspan(t+1))) < toler)
                
                tic

                tspan_bk = tspan(1:t+1);
                sol = sol(:,1:t+1);

                time=time+toc;
                break
            end
            
            tic

            tspan_bk = tspan;
            
            time=time+toc;
        end
    
    end
    
end