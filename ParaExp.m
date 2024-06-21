function [X_para,time] = ParaExp(p,tspan,M,G,x0)

    Type_2_expm = true;

    time = 0;

    mx              = length(x0);
    nt              = length(tspan);
    size_int        = fix(nt/p);

    X_para      = zeros(mx,nt);
    X_para(:,1) = x0;

    ini_zero        = zeros(mx,1);

    x0p            = zeros(mx,p);
    x0p(:,1)       = x0;

    for iter = 1:p
    
    %----- Interval definition
        
        deb       = 1*(iter==1) + (iter-1)*size_int*(iter > 1);
        fin       = iter*size_int;
        if iter == p
            fin = length(tspan);
        end
        time_interval = tspan(deb:fin);
        
    %----- Type 1
                
        [~,X_T1_p] = odeRK4_inhom_ufunc(M,G,time_interval,ini_zero);        
        X_para(:,deb:fin) = X_para(:,deb:fin) + X_T1_p;
        
        if iter < p
            x0p(:,iter+1) = X_T1_p(:,end);
        end

    end

    for iter = 1:p
    
    %----- Interval definition
        
        deb       = 1*(iter==1) + (iter-1)*size_int*(iter > 1);
        if Type_2_expm
            time_interval_till_end_from_zero = tspan(deb:end) - tspan(deb);
        else
            time_interval_till_end = tspan(deb:end);
        end

        
    %----- Type 2

        if Type_2_expm
            [X_T2_p] = expm_prop(M,time_interval_till_end_from_zero,x0p(:,iter));
        else
            [~,X_T2_p] = odeRK4_inhom_ufunc(M,@(t) zeros(mx,1),time_interval_till_end_from_zero,x0p(:,iter));
        end
                
        X_para(:,deb+1:end) = X_para(:,deb+1:end) + X_T2_p(:,2:end);
                
    end

end