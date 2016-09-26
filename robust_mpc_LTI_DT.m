function robust_mpc_LTI_DT
% example_8_31
% This example implements Artstein's circles. The boundary for the state
% solution can be shifted by setting the variable bound. Moreover, both the
% continuous time and the discrete time version are implemented and can be
% switched by supplying "system_ct" instead of "system_dt" to the NMPC routine.
% Note that the variable "type" needs to be reset as well.
    addpath('./nmpcroutine');
    close all;

    global psys real_psys gamma beta;
    
    % System definition
    T    = 1;
    csys = passive_ss([1 1 1]);          % SOP and ZSD system
    real_csys = passive_ss([1 1 2.9]);   % SOP and ZSD system
    psys = c2d_passive_system(csys, T);
    real_psys = c2d_passive_system(real_csys, T);
    
    % MPC parameters
    mpciterations = 20;
    N             = 4;
    gamma         = 1;
    beta          = 1;
    t0            = 0.0;
    
    % Initial conditions
    x0      = [2.4, -1.8];
    z0      = 0;
    x_ext0  = [x0, z0];
    %u_      = 1;
    u0      = ones(1,N);
    
    % Tolerances and parameters
    tol_opt       = 1e-8;
    opt_option    = 0;
    iprint        = 5;
    type          = 'difference equation';
    atol_ode_real = 1e-12;
    rtol_ode_real = 1e-12;
    atol_ode_sim  = 1e-4;
    rtol_ode_sim  = 1e-4;
    
%     [beta, gamma] = find_params( x_ext0, u0, N, t0, T, ...
%              tol_opt, opt_option, type, ...
%              atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
%              iprint)

    nmpc(@runningcosts, @terminalcosts, @constraints, ...
         @terminalconstraints, @controllinearconstraints, ...
         @numlinearconstraints, @system_dt, @realsystem_dt, ...
         mpciterations, N, T, t0, x_ext0, u0, ...
         tol_opt, opt_option, ...
         type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
         iprint, @printHeader, @printClosedloopData, @plotTrajectories);

    rmpath('./nmpcroutine');
end


function [ beta, gamma ] = find_params( x_ext0, u0, N, t0, T, ...
             tol_opt, opt_option, type, ...
             atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
             iprint)
         
    p = size(u0,2);
    new_u0 = [ones(1,p); ones(1,p); u0];
    
    new_sys = @(t, x_ext, u, T) system_dt(t, x_ext, u(3:end), T);
    new_running_costs = @(t, x_ext, u) runningcosts(t, x_ext, u(3:end));
    
    [~ ,~ ,u] = nmpc(new_running_costs, @terminalcosts, @constraints, ...
         @terminalconstraints, @params_controllinearconstraints, ...
         @numlinearconstraints, new_sys, new_sys, ...
         1, N, T, t0, x_ext0, new_u0, ...
         tol_opt, opt_option, ...
         type, atol_ode_real, rtol_ode_real, atol_ode_sim, rtol_ode_sim, ...
         iprint)
     
     beta  = u(1);
     gamma = u(2);
end

function [A, b, Aeq, beq, lb, ub] = params_controllinearconstraints(t, x, u)

    A   = -eye(2,length(u(:)));
    b   = zeros(2,1);
    Aeq = [];
    beq = [];
    lb  = -1;
    ub  =  1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of the NMPC functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cost = runningcosts(t, x_ext, u)
    % Extract info from extended state
    x = x_ext(1:length(x_ext)-1)';
    u = u';
    
    cost = 2*(x')*x + u*u;
end

function cost = terminalcosts(t, x_ext)
    % Extract info from extended state
    
    x = x_ext(1:end-1)';
    
    cost = 2*(x')*x;
end

function [c,ceq] = constraints(t, x, u)
    c   = [];
    ceq = [];
end

function [c,ceq] = terminalconstraints(t, x_ext)
    global beta;
    % Extract info from extended state
    x = x_ext(1:length(x_ext)-1)';
    z = x_ext(length(x_ext))';
    
    c(1) = z - beta;
    ceq = x';
end

function [Nconstr, Nconstr_eq] = numlinearconstraints(constraints_fun, t, x, u)
    [~, b, ~, beq, ~, ~] = constraints_fun(t, x, u);
    Nconstr = size(b,1);
    Nconstr_eq = size(beq,1);
end

function [A, b, Aeq, beq, lb, ub] = controllinearconstraints(t, x, u)
    A   = [];
    b   = [];
    Aeq = [];
    beq = [];
    lb  = -1;
    ub  =  1;
end

function x_ext_new = system_dt(t, x_ext, u, T)
    global psys gamma;
    % Extract info from extended state
    x = x_ext(1:length(x_ext)-1)';
    z = x_ext(length(x_ext))';
    u = u';
    % State equations
    x_new = psys.A*x + psys.B*u;
    y = psys.C*x + psys.D*u;
    z_new = -u'*y + gamma*(u')*u + z;   % Ensuring controller SOP
    % Next extended state
    x_ext_new = [x_new', z_new'];
end

function x_ext_new = realsystem_dt(t, x_ext, u, T)
    global real_psys gamma;
    
    psys = real_psys;
    % Extract info from extended state
    x = x_ext(1:length(x_ext)-1)';
    z = x_ext(length(x_ext))';
    u = u';
    % State equations
    x_new = psys.A*x + psys.B*u;
    y = psys.C*x + psys.D*u;
    z_new = -u'*y + gamma*(u')*u + z;   % Ensuring controller SOP
    % Next extended state
    x_ext_new = [x_new', z_new'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of output format
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printHeader()
    fprintf('   k  |      u(k)        x(1)        x(2)     Time\n');
    fprintf('--------------------------------------------------\n');
end

function printClosedloopData(mpciter, u, x, t_Elapsed)
    fprintf(' %3d  | %+11.6f %+11.6f %+11.6f  %+6.3f', mpciter, ...
            u(1,1), x(1), x(2), t_Elapsed);
end

function plotTrajectories(dynamic, system, real_system, T, t0, x0, u, ...
                          atol_ode, rtol_ode, type)
    [~, t_intermediate, x_intermediate] = dynamic(system, T, t0, ...
                                          x0, u, atol_ode, rtol_ode, type);
                                      
    [~, real_t_intermediate, real_x_intermediate] = dynamic(real_system, T, t0, ...
                                          x0, u, atol_ode, rtol_ode, type);
    figure(1);
        title('x_1/x_2 closed loop trajectory');
        xlabel('x_1');
        ylabel('x_2');
        grid on;
        hold on;

        plot(x_intermediate(:,1),x_intermediate(:,2),'-og', ...
             'MarkerFaceColor','g');
        plot(real_x_intermediate(:,1),real_x_intermediate(:,2),'-or', ...
             'MarkerFaceColor','r');
        axis([-3 3 -3 3]);
        axis square;

    figure(2);
        title('x_1 and x_2 closed loop trajectory');
        xlabel('n');
        ylabel('x_1(n), x_2(n)');
        grid on;
        hold on;
        plot(t_intermediate,x_intermediate(:,1),'-og');
        plot(t_intermediate,x_intermediate(:,2),'-og');
        
        plot(real_t_intermediate,real_x_intermediate(:,1),'-or');
        plot(real_t_intermediate,real_x_intermediate(:,2),'-or');
        axis([0 30 -3 3]);
        axis square;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
