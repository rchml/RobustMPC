function results = gurobi_qcqp( Q, c, alpha, A, b, lsense, lb, ub, Qc, q, beta, init, printInfo )

% gurobi_qcqp   solves a quadratic constrained quadratic problem using
%               the gurobi solver.
% 
%  QCQP model:
%  minimize
%      x'Qx + cx + alpha
%  subject to
%      Ax 'lsense(>,<,=)' b
%      x'Qcx + q'x   <=   beta
%      lb   <=   x   <=   ub
%      x0(i)          =   init(i)   \forall i
%      x(i) \in R                   \forall i

    clear model;
    
    % Objective function
    model.Q      = sparse(Q);
    model.obj    = c;
    model.objcon = alpha;   
    % Linear constraints
    model.A      = sparse(A);
    model.rhs    = b';
    model.sense  = lsense;
    % Bound constraints
    if ~isempty(lb)
        model.lb     = lb;
    end
    if ~isempty(ub)
        model.ub     = ub;
    end
    
    % Quadratic constraints
    for i = 1:length(Qc)
        model.quadcon(i).Qc  = sparse(Qc{i});
        model.quadcon(i).q   = q{i};
        model.quadcon(i).rhs = beta{i};
    end
    % Initial values
    model.start = init';
    
    %gurobi_write(model, 'qcqp.lp');
    if printInfo
        params.OutputFlag   = 1;
    else
        params.OutputFlag   = 0;
    end
    
    params.ScaleFlag        = 1;
    params.NumericFocus     = 3; 
    params.BarHomogeneous   = 1;
    params.DualReductions   = 0;
    params.quad             = 1;
    
    results = gurobi(model, params);
    
end
