
function [Qz, qz, beta_] = constraints( sys, N, K, beta, gamma_, x0_, x_ext0 )

n = size(sys.B,1);
p = size(sys.B,2);

x0 = sym('x0', [n,1]);
u0 = sym('u0', [p,1]);
u  = sym('u', [p,N-1]);
assume(x0,'real')
assume(u0,'real');
assume(u,'real');

A  = sym('A', [n,n]);
B  = sym('B', [n,p]);
C  = sym('C', [p,n]);
D  = sym('D', [p,p]);

gamma = sym('gamma');
z0    = sym('z0');
    
l = z(K);

l = subs(l,A(:), sys.A(:));
l = subs(l,B(:), sys.B(:));
l = subs(l,C(:), sys.C(:));
l = subs(l,D(:), sys.D(:));
l = subs(l,gamma, gamma_);
l = subs(l, x0, x0_);
u_t = u';
u_var = [u0(:); u_t(:)];

% Calculate Qz
h = 1/2*hessian(l, u_var);
Qz = double(h);

% Calculate qz
g = gradient(l, u_var);
g = subs(g, u_var, zeros(N,1));
qz = double(g);

% Calculate beta_
const = subs(l, u_var, zeros(N,1));
const = subs(const, z0, x_ext0);
beta_ = beta - const;
    

    function out = z(k)

        if k == 1
            out = z0 - u0'*y(0) + gamma*(u0')*u0;
            return
        end
        out = z(k-1) - transpose(u(:,k-1))*y(k-1) + gamma*transpose(u(:,k-1))*u(:,k-1);
        
    end

    function out = y(k)

        if k == 0
            out = C*x0 + D*u0;
            return
        end
        out = C*x(k) + D*u(:,k);
    end

    function out = x(k)

        if k == 1
            out = A*x0 + B*u0;
            return
        end
        out = A*x(k-1) + B*u(:,k-1);
        
    end

end
