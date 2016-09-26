function sys = passive_ss( vect )
% passive_ss returns a passive LTI system ss.
%   vect: should contain a vector of non-negative numbers representing
%             the inductance, capacitance and resistance values of a
%             parallel circuit of series-connected bipoles. 
%             The ss obj returned is a minimal realization of the
%             current/voltage transfer function of the circuit. 
%
%             Input are expected in this format: 
%             [L1, C1, R1, L2, C2, R2,...]
%             Then vect must be s.t. length(vect) = 3*k, with k = 1,2,... 

    assert(mod(length(vect), 3) == 0, ['Error: vect should be s.t.'...
                        'length(vect) = 3*k, with k = 1,2,...']); 
    n = length(vect)/3;
    tmp_tf = tf(0,1);
    for i = 1:n
        l = vect(3*i-2);
        c = vect(3*i-1);
        r = vect(3*i);
        tmp_tf = tmp_tf + tf([c 0], [l*c r*c 1]);
    end
    sys = ss(tmp_tf);
end

