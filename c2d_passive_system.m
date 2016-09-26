function [ dsys ] = c2d_passive_system( csys, T )
% c2d_passive_system( csys, T ) 
% Returns a discrete passive system from a continuous time passive one.
% csys:     a ss continuous time passive system with D = 0 and
%           dim(Ker(A)) = 0
% T:        sample time of the discretization

assert (rank(csys.A) == size(csys.A,1))

    n = size(csys.A, 1);
    invAxB = csys.A \ csys.B;
    expAT = expm(csys.A*T);
    
    Ad = expAT;
    Bd = (expAT - eye(n,n)) * invAxB;
    Cd = csys.C * (csys.A \ (expAT - eye(n,n)));
    Dd = csys.C * ((expAT - eye(n,n)) * (csys.A*csys.A) \ csys.B - invAxB * T);

    dsys = ss(Ad, Bd, Cd, Dd, T);

end

