function u = LyapunovController(roe, roe_des, oe_c, params)
    % Compute Lyapunov control input (roes are all mean)
    % Returns u
    
    % unpack parameters: [N, k, u_lowerbound, u_upperbound, dlambda_thresh, dlambda_dot] = params;
    N = params.N;
    k = params.k;
    u_lowerbound = params.u_lowerbound;
    u_upperbound = params.u_upperbound;
    dlambda_thresh = params.dlambda_thresh;
    dlambda_dot = params.dlambda_dot;
    
    Delta_roe = roe - roe_des;
    if abs(Delta_roe(2)) >= dlambda_thresh
        % impose a desired da to cause a desired drift in dlambda
        mu_E = 3.986e5;
        a_c = oe_c(1);
        n_c = sqrt(mu_E / a_c^3);
        da_des = -2/3 * sign(-Delta_roe(2)) * dlambda_dot / n_c;
        Delta_roe(1) = roe(1) - da_des;
    end
    Delta_roe_reducedmodel = [Delta_roe(1); Delta_roe(3:6)];
    roe_d_mean_reducedmodel = [roe(1); roe(3:6)];
    
    A = A_Control_ReducedModel(oe_c(:));
    B = B_Control_ReducedModel(oe_c(:));
    P = P_Control_ReducedModel(oe_c(:), Delta_roe_reducedmodel, N, k);
    u = - pinv(B) * (A * roe_d_mean_reducedmodel + P * Delta_roe_reducedmodel);
    
    % check lower bound for control effort
    if abs(u(1)) < u_lowerbound
        u(1) = 0;
    end
    if abs(u(2)) < u_lowerbound
        u(2) = 0;
    end
    % check upper bound for control effort
    if abs(u(1)) > u_upperbound
        u(1) = sign(u(1)) * u_upperbound;
    end
    if abs(u(2)) > u_upperbound
        u(2) = sign(u(2)) * u_upperbound;
    end
    
    u = [0; u];
end

