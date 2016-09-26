
function [x_tot, code, beta, failureCount] = MPC_robust_LTI_DT( sys, ...
                                real_sys, x0, T, Q, R, P, N, mpc_iter, ...
                                betaPolicy, beta, P_beta, gamma_, z_fun, c_path, hard_reset )
                            
    % code = 0  : OK
    % code = -1 : couldn't find beta
    % code = -2 : couldn't find u
    % code = -3 : Invalid z_fun
    % code = -4 : huge value of state vars

    nominal_sim     = 0;
    relax_mod       = 0;
    forceSymCalc    = 0;
    doPrints        = 1;
    doPlot          = 1;
    betaMonotone    = 0;
    
    p = size(sys.B, 2);
    n = size(sys.B, 1);

    % Initial state (controller)
    x_ext0 = 0;
       
    % Initial control vector
    z0 = zeros(1,N)';
    
    % Solve the quadratic constrained optimization problem
    quadratic_con = 1;
    
    % Constraints parameters
    % Selecting gamma in order to have Qc PSD and the selected condition satisfied
    if ~isempty(gamma_)
        if strcmp(z_fun, 'SOP')
            gamma_tmp = max(real(eig(-sys.D)), 1/gamma_);   % 1/gamma : upper bound to the controller L2-gain                         
        elseif strcmp(z_fun, 'L2')
            gamma_tmp = gamma_;                             % gamma   : upper bound to the controller L2-gain
        end
    else
        gamma_tmp = max(real(eig(-sys.D)), 0);
    end
        
    % Define beta_(x0) = x0' P x0
    P_beta = get_lyap_matrix(sys);
    beta_ = @(x0, z0) (x0')*P_beta*x0 + z0;
    
    % Load or create problem struct
    mpc_prob = get_problem_data(sys, Q, P, R, N, gamma_tmp, quadratic_con, z_fun, forceSymCalc);
    
    % Get calculated gamma value
    gamma = mpc_prob.gamma;
    
    % Initialize failureCount, num of failure in finding an opt control,
    % which is replaced by an almost feasible control when possible
    failureCount = 0;
    
    % State trajectory vector
    x_tot = zeros(n, mpc_iter+1);
    x_tot(:,1) = x0;
    
    % Expected state trajectory vector
    x_exp_tot = zeros(n, mpc_iter);
    
    % Selected control vector
    z_tot = zeros(p, mpc_iter);
    
    % Selected beta vector
    beta_vect = zeros(1, mpc_iter);
    
    if nominal_sim
        real_sys = sys;
    end

    for i = 1:mpc_iter
        
        % Exit condition: huge value of state variables
        if norm([x0; x_ext0], Inf) > 1e15
            disp('WARNING: huge value of the state variables, simulation STOPPED\n');
            % Print data and plots
            save_data( n, p, T, N, z_tot(:,1:i), x_tot(:,1:i), ...
                             x_exp_tot(:,1:i), beta_vect(1:i), ...
                                        doPrints, doPlot, i, c_path );
            x_tot = [];
            code  = -4;
            fprintf('Iteration: %d\n', i);
            return
        end
        
        % Step (1) Obtain initial state (feedback)
        %x0 = measure_initial_state(real_sys, x0);
        
%         % Do not re-compute the same symbolical calculation at each execution
%         if hard_reset
%             reset = (i == 1);
%         else
%             reset = 0;
%         end
                
        % Step (1.5) Calculate new beta
        
        switch betaPolicy
            case 0
                beta = beta_(x0, x_ext0);

            case 1
                z0_beta = [0; z0];
                
                beta_lb = 0;
                if betaMonotone
                    beta_lb = beta;
                end
                
                [beta, z0]  = get_beta_feas(mpc_prob, x0, x_ext0, z0_beta, beta_lb);
                beta = beta+1;                                                         %%%%%%%
            case 2
                z0_beta = [0; z0];
                
                beta_lb = 0;
                if betaMonotone
                    beta_lb = beta;
                end
                
                [beta0, z0] = get_beta_feas(mpc_prob, x0, x_ext0, z0_beta, beta_lb);
                beta = max([beta_(x0, x_ext0), beta0]);
                
            otherwise
                disp('ERROR: betaPolicy INVALID\n');
                code = -1;
                return
        end
        
       
        if isempty(beta) || (betaPolicy == 2 && isempty(beta0))
            % Print data and plots
            save_data( n, p, T, N, z_tot(:,1:i), x_tot(:,1:i), ...
                             x_exp_tot(:,1:i), beta_vect(1:i), ...
                                        doPrints, doPlot, i, c_path );
            x_tot = [];
            code  = -1;
            fprintf('Iteration: %d\n', i);
            return
        end
        
        % Update total beta vector
        beta_vect(i) = beta;
        
        % Step (2) Solve the optimal control problem
        [u, failureCount] = mpc_robust_control(mpc_prob, x0, x_ext0, z0, ...
                                        beta, betaPolicy, failureCount );
        
        if isempty(u)
            
            % Try once more removing quadratic constraints
            if relax_mod
%                 quadratic_con = 0;
                [u, failureCount] = mpc_robust_control(mpc_prob, x0, x_ext0, z0, ...
                                        beta, betaPolicy, failureCount );
            end
            
            if isempty(u)
                % Print data and plots
                save_data( n, p, T, N, z_tot(:,1:i), x_tot(:,1:i), ...
                                x_exp_tot(:,1:i), beta_vect(1:i), ...
                                            doPrints, doPlot, i, c_path );
                x_tot = [];
                code  = -2;
                fprintf('Iteration: %d\n', i);
                return
            end
        end
        
        % Step (3) Apply control to the system
        x_next      = real_sys.A*x0 + real_sys.B*u;     % Sys state
        y           = real_sys.C*x0 + real_sys.D*u;     % Sys output
        
        % Controller state
        if strcmp(z_fun, 'SOP')
            x_ext_next  = (u')*y + gamma*(u')*u + x_ext0;
        elseif strcmp(z_fun, 'L2')
            x_ext_next  = (u')*u - gamma*(y')*y + x_ext0;
        else
            code = -3;
            return;
        end
        
        % Update total state vector
        x_tot(:, i+1) = x_next;
        
        % Update total control vector
        z_tot(:, i) = u;
        
        % Expected state
        x_exp_next = sys.A*x0 + sys.B*u;
        
        % Update total expected state vector
        x_exp_tot(:, i) = x_exp_next;
                
        % Update initial state
        x0     = x_next;
        
        % Update controller initial state
        x_ext0 = x_ext_next;
    end
    
    code = 0;
       
    % Print data and plots
    save_data( n, p, T, N, z_tot, x_tot, x_exp_tot, beta_vect, ...
                                    doPrints, doPlot, mpc_iter, c_path );
        
end

function mpc_struct = get_problem_data( sys, Q, P, R, N, gamma, quadratic_con, z_fun, forceSymCalc )

    if ~forceSymCalc && exist(fullfile(cd, 'mpc_prob.mat'), 'file')
        
        mpc_prob = load('mpc_prob.mat');

        if ~isempty(mpc_prob) && ...
                isequal(mpc_prob.sys, sys) && ...
                isequal(mpc_prob.Q, Q) && ...
                isequal(mpc_prob.P, P) && ...
                isequal(mpc_prob.R, R) && ...
                isequal(mpc_prob.N, N) && ...
                isequal(mpc_prob.z_fun, z_fun)

            mpc_struct = mpc_prob;
        else
            mpc_struct = qcqp_mpc_data( sys, Q, P, R, N, gamma, quadratic_con, z_fun);
        end
    else
        mpc_struct = qcqp_mpc_data( sys, Q, P, R, N, gamma, quadratic_con, z_fun);
    end
    
end

function save_data( n, p, T, N, z_tot, x_tot, x_exp_tot, beta_vect, ...
                        doPrints, doPlot, n_iter, c_path )

    sim_ended = size(x_tot,2) > n_iter;
    
    % Plots and prints
    if doPrints
        s = printHeader(n,p, [], N);
        for i = 1:n_iter
            s = strcat(s, printClosedLoopData(i,z_tot(:,i),x_tot(:,i), beta_vect(i)));
        end
        if sim_ended
            s = strcat(s, printFinalState(x_tot(:,n_iter + 1)));
        end
        
        % Write output.txt
        fid = fopen(strcat(c_path, sprintf('/output.txt')), 'wt');
        fprintf(fid, s);
        fclose(fid);
        
        % Save simulation variables
        data.x_tot      = x_tot;
        data.beta_vect  = beta_vect;
        data.u          = z_tot;
        save(strcat(c_path, sprintf('/simdata.mat')), '-struct', 'data');

        if size(x_tot, 2) == 1
            return
        end
        
        if doPlot
            if n == 2
                fh1 = plot2DStateTrajectories(x_exp_tot, x_tot);
                saveas(fh1,strcat(c_path, sprintf('/x_wpred.png')));
                close(fh1);
            end
            if sim_ended
                fh2 = plotTrajectories(x_tot, T*(1:n_iter));
            else
                fh2 = plotTrajectories(x_tot, T*(1:n_iter-1));
            end
            
            % Plotting Beta
            fh3 = figure();
            plot(T*(1:n_iter), beta_vect);
            grid on;
            title('Choosen Beta');
            xlabel('time');
            ylabel('Beta');
            
            saveas(fh2,strcat(c_path, sprintf('/xdata.png')));
            savefig(fh2, strcat(c_path, sprintf('/xdata')));
            saveas(fh3,strcat(c_path, sprintf('/Betadata.png')));
            
            close(fh2);
            close(fh3);
        end
    end
end


function [u, failureCount] = mpc_robust_control( mpc_prob, x0, x_ext0, z0, ...
                                        beta, betaPolicy, failureCount )

    p = size(mpc_prob.sys.B, 2);
    n = size(mpc_prob.sys.B, 1);
    
    printInfo   = 0;
    z0_empty    = 0;
    
    if isempty(z0)
        z0_empty    = 1;
        z0          = zeros(mpc_prob.N*p, 1);
    end
    
    % Do not re-compute the same symbolical calculation at each execution
    [Q_o, c_o, alpha_o, A, b, lsense, lb, ub, Qc, q, beta_c]...
            = get_qcqp_mpc_matrices( mpc_prob, x0, x_ext0, beta );
    init = z0;
        
%     % Do not re-compute the same symbolical calculation at each execution
%     [Q_o, c_o, alpha_o, A, b, lsense, lb, ub, Qc, q, beta_c]...
%             = qcqp_mpc_2( sys, x0, x_ext0, Q, P, R, N, gamma, beta, quadratic_con, reset );
%     
%     init = [z0; zeros((N+1)*n + N*p, 1)];
        
    results = gurobi_qcqp(Q_o, c_o, alpha_o,  ...
                          A, b, lsense, ...
                          lb, ub,       ...
                          Qc, q, beta_c,  ...
                          init, printInfo);
                      
    if exist('results', 'var') 
        if isfield(results, 'x')
            u = results.x(1:p);
            return
        end
    end
    
    if betaPolicy == 1 && ~z0_empty
        u = z0(1:p);
        failureCount = failureCount + 1;
        return
    end
    
    fprintf('\nCONTROL: Error code: %s', results.status);
    u = [];
           
end

% Compute the P matrix of the energy function the simple negative output 
% feedback controller
function P = get_lyap_matrix( sys )

    A = sys.A;
    B = sys.B;
    C = sys.C;
    D = sys.D;
    I = eye(size(D,1));
    
    C_ = (I + D)\C;
    A_ = A - B*C_;
    Q_ = (C_')*C_;
    
    P = dlyap(A_', Q_);    
end

function [beta, z] = get_beta_feas( mpc_prob, x0, x_ext0, z0, beta_lb )

    [Q_o, c_o, alpha_o, A, b, lsense, lb, ub, Qc, q, beta_c]...
            = qcqp_mpc_beta(mpc_prob, x0, x_ext0, beta_lb);
    
    printInfo = 0;
    
    p = size(mpc_prob.sys.B, 2);
    n = size(mpc_prob.sys.B, 1);
    
    init = z0;
%     init = [z0; zeros((N+1)*n + N*p, 1)];
    
    results = gurobi_qcqp(Q_o, c_o, alpha_o,  ...
                          A, b, lsense, ...
                          lb, ub,       ...
                          Qc, q, beta_c,  ...
                          init, printInfo);
   
    if exist('results', 'var') 
        if isfield(results, 'x') 
            beta    = results.x(1);                                            % Constraint check !!!!!!!!!!!!!!!!!!
            z       = results.x(2:end);
            return
        end 
    end
    
    fprintf('\nBETA: Error code: %s\n', results.status);
    beta    = [];
    z       = [];
end
    
function str = printHeader(n, p, beta, N)
    
    if ~isempty(beta)
        str = sprintf('\nChosen beta = %5.4f \nPrediction window size N = %d\n', beta, N);
    else
        str = sprintf('\nChosen beta = variable \nPrediction window size N = %d\n', N);
    end
        
%     pause;
    
    u = @(k) sprintf('    u(%d)    ',k);
    x = @(k) sprintf('    x(%d)    ',k);
    
    s = '\n     k |  ';
    for i = 1:p
        s = strcat(s, u(i));
    end
    
    for i = 1:n
        s = strcat(s, x(i));
    end
    
    if isempty(beta)
        s = strcat(s, '    beta    ');
    end
    
    s = strcat(s, '\n----------------------------------------------------------------------------------\n');
        
    str = strcat(str, s);
%     fprintf(s);
end

function s = printFinalState(x)
    
    x_s = @(k) sprintf('  x(%d) = ',k);
    
    n = size(x,1);
    s = '\nLast state: ';

    for i = 1:n
        s = strcat(s, x_s(i));
        s = strcat(s, sprintf(' %5.4f   ', x(i)));
    end
    
    s = strcat(s, '\n----------------------------------------------------------------------------------\n');
    
%     fprintf(s);
end

function printBetaVector(beta_vect)

    b = size(beta_vect,2);
    if ~isempty(beta_vect)
        for i = 1:b
            s = sprintf(' beta(%d) =   %2.8f', i, beta_vect(i));
            fprintf([s, '\n']);
        end
    end
end

function s = printClosedLoopData(mpciter, u, x, beta)
    
    p = size(u,1);
    n = size(x,1);
    
    s = sprintf('%6d   ', mpciter);
    for i = 1:p
        s = strcat(s, [sprintf('   %5.4f  ', u(i)),' ']);
    end
    for i = 1:n
        s = strcat(s, sprintf('   %5.4f   ', x(i)));
    end
    
    if ~isempty(beta)
            s = strcat(s, sprintf('   %2.8f   ', beta));
    end
    
    s = strcat(s, '\n');
%     fprintf([s, '\n']);
end
    


function fh = plot2DStateTrajectories(x_exp, x_real)

    n = size(x_real,1);
    
    assert(n == 2);
    
    fh = figure(1);
        title('x_1/x_2 closed loop trajectory');
        xlabel('x_1');
        ylabel('x_2');
        grid on;
        hold on;

        for i = 1:size(x_exp,2)
            p_1 = plot([x_real(1,i),x_exp(1,i)], [x_real(2,i),x_exp(2,i)],'-og', ...
                 'MarkerFaceColor','g');
        end
        
        p_2 = plot(x_real(1,:),x_real(2,:),'-or', ...
             'MarkerFaceColor','r');
         
        x_max_1 = max(max([x_real(1,:), x_exp(1,:)],[],2),[],1);
        x_min_1 = min(min([x_real(1,:), x_exp(1,:)],[],2),[],1);
        x_max_2 = max(max([x_real(2,:), x_exp(2,:)],[],2),[],1);
        x_min_2 = min(min([x_real(2,:), x_exp(2,:)],[],2),[],1);
        const = 0.1;
        axis([x_min_1-const x_max_1+const x_min_2-const x_max_2+const]);
        axis square;
        
        legend([p_1 p_2], 'Predicted trajectory', 'Actual trajectory');
        legend('show')

end
    

function fh = plotTrajectories(x, t_vect)

    n = size(x,1);
    
    fh = figure();
        title('x closed loop trajectory');
        xlabel('time');
        ylabel('x');
        grid on;
        hold on;
        
        % In order to write the legend
        x_ = @(k) sprintf('x_%d', k);
        
        x_cell = cell(n,1);

        for i = 1:n
            plot([0, t_vect],x(i,:),'-o');
            x_cell(i) = {x_(i)};
        end
        
        axis([0 t_vect(end)+1 min(x(:),[],1) max(x(:),[],1)]);
        axis square;
        
        legend(x_cell{:});
        legend('show')

end
    
    