
% LQ and passive MPC comparison
function LQ_passiveMPC_comparison            
            
    clear all
        
    % Perturbation definition
%     num = 10;
%     k   = 8;
%     x   = k:num + k-1;
%     y   = k:num + k-1;
%     ax  = .4;
%     ay  = .3;
        
    % IQC choice ( 0: passivity, 1: small gain)
    IQC = 0;
            
    % Policy for Beta e Gamma parameters
    betaPol = 1;
    gamPol  = 1;
    
    % System definition
%     T   = 1;
%     R1  = .1;
%     L1  = .8;
%     L2  = .5;
%     R2  = .2;
%     xn  = .5;
%     yn  = 3;
    
    T   = 1;
    R1  = .1;
    L2  = .5;
    C2  = 3;
    R2  = .2;
    xn  = .8;
    yn  = .2;
    
    num = 1;
    k   = 1;
    x   = k:num + k-1;
    y   = k:num + k-1;
    ax  = .2;% - xn; %1.2;
    ay  = .7;% - yn; %.6;

                                              
%     csys = passive_ss([1,1,1]);
    spaceDim = 4;
    csys = passive_ss([xn, yn, R1, L2, C2, R2]);             % SOP and ZSD system
%     csys = passive_ss([L1, xn, R1, L2, yn, R2]);             % SOP and ZSD system
    sys = c2d_passive_system(csys, T);
    
    % Create perturbed sys set
    [p_sys_set, p_ind_set] = set_of_perturbed_sys(T, xn, yn, R1, L2, C2, R2, num, x, y, ax, ay);
%     [p_sys_set, p_ind_set] = set_of_perturbed_sys(T, xn, yn, R1, L1, L2, R2, num, x, y, ax, ay);
    
    p = size(sys.B, 2);
    n = size(sys.B, 1);
    
    % Optimal control weights
    Q = 2.3*eye(n);
    R = 1.2*eye(p);
        
    % Window size
    N  = 10; 
    
    x0 = .1*ones(n,1);
    mpc_iter = 1000;
    threshold = .05;
    interval = 3;
    
    if IQC == 0
        z_fun  = 'SOP';
    else 
        z_fun  = 'L2';
        gamPol = 1;
    end
    
    xmin = ax*x(1)+xn;
    xmax = ax*x(num)+xn;
    ymin = ay*y(1)+yn;
    ymax = ay*y(num)+yn;
    
    % Main path definition
    m_path = sprintf('DATA_SIM/exp_xR%d_N%d_C1=%0.1f[%0.1f,%0.1f]_C2=%0.1f[%0.1f,%0.1f]_BetaPol%d_GamPol%d_%s_beta+10', ...
            spaceDim, N, xn, xmin, xmax, yn, ymin, ymax, betaPol, gamPol, z_fun);
        
        tic
    % Calculate stability logical matrix for the lq closed loops
    mpc_stability_matrix = check_mpc_robustness(sys, p_sys_set, p_ind_set, ...
                                      x0, T, Q, R, N, ...
                                      mpc_iter, gamPol, betaPol, z_fun, ...
                                      interval, threshold, m_path);
                                  
        toc
                                           
    % Calculate stability logical matrix for the lq closed loops
    lq_stability_matrix = check_lq_robustness(sys, p_sys_set, Q, R, N);
    
    % Draw stability plot for the mpc closed loops
    fh = figure('units', 'normalized', 'position', [.1 .1 .7 .7]);
    n_plot = 1;
    draw_stability_plot(mpc_stability_matrix, ax, ay, xn, yn, x, y, n_plot);
        
    % Draw stability plot for the lq closed loops
    n_plot = 2;
    draw_stability_plot(lq_stability_matrix, ax, ay, xn, yn, x, y, n_plot);
    
    saveas(fh,strcat(m_path, sprintf('/stability_plot.png')));
    savefig(fh, strcat(m_path, sprintf('/stability_plot')));
    close(fh);

end

% Compute the P matrix of the energy function of the simple negative output 
% feedback controller
function P = get_lyap_matrix( sys )

    k = 1;
    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;
    I = eye(size(D,1));
    
    C_ = (I + D*k)\C;
    A_ = A - B*k*C_;
    Q_ = (C_')*C_;
    
    P = dlyap(A_', Q_);    
end

function stability_matrix = check_mpc_robustness( sys, p_sys_set, p_ind_set, ...
                                        x0, T, Q, R, N, ...
                                        mpc_iter, gamPol, betaPol, z_fun, ...
                                        interval, threshold, m_path )

    num = size(p_sys_set,1);
    P = Q;
    
    stability_matrix = zeros(num, num);
    
    % Do not re-compute the same symbolical calculation at each execution
    reset           = 1;
    beta            = 0;
    P_beta = get_lyap_matrix(sys);
    
    delta = 0;
    
    if gamPol == 1 || strcmp(z_fun, 'L2')
        for i = 1:num
            for j = 1:num
                p_sys = p_sys_set{i,j};
                % if gamPol == 1 select gamma as an upper bound of the Hinf norm of the p_sys
                delta = max([delta, norm(p_sys, inf)]);
            end
        end
        gamma = 0.9/delta;
    else 
        gamma = [];
    end
    
    % Store mpc_prob.mat
%     quadratic_con = 1;
%     qcqp_mpc_data( sys, Q, P, R, N, gamma, quadratic_con, z_fun);
    
    parfor i = 1:num
        for j = 1:num
            p_sys = p_sys_set{i,j};
                        
            p_ind = p_ind_set{i,j};
            c_path = strcat(m_path, sprintf('/C1=%0.1f_C2=%0.1f', p_ind.x, p_ind.y));
            mkdir(c_path);
            
%             [x_tot, code, beta, failureCount] = MPC_robust_LTI_DT(sys, ...
%                     p_sys, x0, T, Q, R, P, N, mpc_iter, betaPol, beta,
%                     P_beta, gamma, c_path, reset); % COMMENTED FOR parfor
            
            [x_tot, code, ~, failureCount] = MPC_robust_LTI_DT(sys, ...
                    p_sys, x0, T, Q, R, P, N, mpc_iter, betaPol, beta, P_beta, gamma, z_fun, c_path, reset);
                
            if failureCount ~= 0
                fprintf('\nFailure count: %d\n', failureCount);
            end
            
            if (code == 0) && ...
                        (check_stability(x_tot, interval, threshold) == 1)
                stability_matrix(i,j) = 1;
            end
            
            if code == -3
                fprintf('Invalid z_fun\n');
            end
            
            fprintf('.');
        end
        fprintf('\n');
    end
end

function out = check_stability( x, interval, threshold )

    if isempty(find(abs(x(:, end-interval+1:end)) >= threshold,1))
        out = 1;
    else
        out = 0;
    end
end

function [p_sys_set, p_ind] = set_of_perturbed_sys( T, xn, yn, R1, L1, L2, R2, num, x, y, ax, ay )

    p_sys_set = cell(num,num);
    p_ind = cell(num,num);
    
    for i = 1:num
        for j = 1:num
            xp = ax*x(i)+xn;
            yp = ay*y(j)+yn;
                        
%             real_csys = passive_ss([1.1 .9 2.9]);
            real_csys = passive_ss([L1, xp, R1, L2, yp, R2]);       % SOP and ZSD system
%             real_csys = passive_ss([xp, yp, R1, L1, L2, R2]);       % SOP and ZSD system
            p_sys_set(i,j) = {c2d_passive_system(real_csys, T)};
            
            S.x = xp;
            S.y = yp;
            p_ind(i,j)     = {S};
        end
    end
end


function draw_stability_plot(stability_matrix, ax, ay, xn, yn, x, y, n_plot)

    [stab_x, stab_y]     = find(stability_matrix > 0);
    [instab_x, instab_y] = find(stability_matrix == 0);
    
    x_green = ax*x(stab_x)+xn;
    y_green = ay*y(stab_y)+yn;
    
    x_red = ax*x(instab_x)+xn;
    y_red = ay*y(instab_y)+yn;

    subplot(1, 2, n_plot);
    grid on;
    hold on;

    scatter(x_green, y_green, [], 'g', 'filled');
    scatter(x_red, y_red, [], 'r', 'filled');
    scatter(xn, yn, [], 'b', 'filled');
    
    if n_plot == 1
        title('MPC with IQC stability plot');
    else
        title('MPC without IQC stability plot');
    end
    
    xlabel('C1 perturbed sys values');
    ylabel('C2 perturbed sys values');
    axis square
    
    %check_system_stability( K, 31, 60, R1, L2, C2, R2, T )
end


function stability_matrix = check_lq_robustness( sys, p_sys_set, Q, R, N )
    
    num = size(p_sys_set,1);
    
    % Calculate linear static control gain
%     K = static_lq_finite_horizon_control( sys, Q, R, N )
    K = static_constrained_lq_finite_horizon_control_SISO( sys, Q, R, N);
    
    stability_matrix = zeros(num, num);
    
    for i = 1:num
        for j = 1:num
            p_sys = p_sys_set{i,j};
            
            if max(abs(eig(p_sys.A + p_sys.B*K))) < 1
                stability_matrix(i,j) = 1;
            end
        end
    end
    
end


function check_system_stability( K, xp, yp, R1, L2, C2, R2, T )

    real_csys = passive_ss([xp, yp, R1, L2, C2, R2]);       % SOP and ZSD system
    real_sys = c2d_passive_system(real_csys, T);
          
    if max(abs(eig(real_sys.A + real_sys.B*K))) < 1
        fprintf('\n xp = %f\n yp = %f\n\nClosed loop is A.S\n', xp, yp);
    else
        fprintf('\n xp = %f\n yp = %f\n\nClosed loop is UNSTABLE\n', xp, yp);
    end
    
end

function K0 = static_constrained_lq_finite_horizon_control_SISO( sys, Q, R, N)

    
    n  = size(sys.B, 1);
    p  = size(sys.B, 2);
    
    K0 = zeros(p, n);
    I  = eye(n);
    
    for i = 1:n
        x0 = I(:,i);

        quadratic_con   = 0;
        printInfo       = 0;
        reset           = 1;

        % Do not re-compute the same symbolical calculation at each execution
        [Q_o, c_o, alpha_o, A, b, lsense, lb, ub, Qc, q, beta_c]...
                = qcqp_mpc( sys, x0, [], Q, Q, R, N, [], [], quadratic_con, reset );

        results = gurobi_qcqp(Q_o, c_o, alpha_o,  ...
                              A, b, lsense, ...
                              lb, ub,       ...
                              Qc, q, beta_c,  ...
                              zeros(N*p,1), printInfo);

        if exist('results', 'var') 
            if isfield(results, 'x')%strcmp(results.status, 'INF_OR_UNBD')
                K0(i) = results.x(1:p);
                continue
            end
        end

        K0 = [];
        return
    end
    
end

function K0 = static_lq_finite_horizon_control( sys, Q, R, N )

    A = sys.A;
    B = sys.B;
    P_N = 10E10*Q;
    
    P0 = dre(P_N, A, B, Q, R, N);
    
    K0 = -(R + B'*P0*B)\(B'*P0*A);

end


function P0 = dre( P_N, A, B, Q, R, N )

    if N == 0
        P0 = P_N;
        return
    end
    
    P_N_prev = A'*P_N*A - (A'*P_N*B) * (inv(R + B'*P_N*B) * (B'*P_N*A)) + Q;
    P0 = dre(P_N_prev, A, B, Q, R, N-1);
end

