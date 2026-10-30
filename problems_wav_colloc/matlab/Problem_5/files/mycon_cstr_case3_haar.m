function [c, ceq] = mycon_cstr_case3_haar(w, m, n_u, Hi, iDWT, z0, tf, u4_fixed)
% Restricoes - Case III (CSTR nao linear), versao Haar (Hi/iDWT).
%
% w = [xd1(1:m), ..., xd8(1:m), u1(1:m), ..., u_nu(1:m)]
%
% Reconstrucao (tempo normalizado tau=[0,1], t=tau*tf):
%   z_i(tau) = z0(i) + (xd_i*inv(Hi))*iDWT     para i=1..8
%   u_i(tau) = amostras diretas (nao integradas)
%
% ceq: residuo da EDO para os 8 estados (fator tf, pois d/dtau=tf*d/dt)
% c:   [] (sem restricoes de desigualdade alem dos limites de u, ja
%       impostos via lb/ub no script principal)

    A = Hi \ iDWT;

    Z = zeros(8, m);
    for i = 1:8
        xdi = w((i-1)*m+1 : i*m);
        Z(i,:) = z0(i) + xdi*A;
    end

    U = zeros(n_u, m);
    for i = 1:n_u
        U(i,:) = w(8*m + (i-1)*m + 1 : 8*m + i*m);
    end

    if n_u == 3
        u1 = U(1,:); u2 = U(2,:); u3 = U(3,:); u4 = u4_fixed*ones(1,m);
    else
        u1 = U(1,:); u2 = U(2,:); u3 = U(3,:); u4 = U(4,:);
    end

    z1=Z(1,:); z2=Z(2,:); z3=Z(3,:); z4=Z(4,:); z5=Z(5,:); z6=Z(6,:); z7=Z(7,:);

    s = u1 + u2 + u4;

    dz1 = u4 - s.*z1 - 17.6*z1.*z2 - 23*z1.*z6.*u3;
    dz2 = u1 - s.*z2 - 17.6*z1.*z2 - 146*z2.*z3;
    dz3 = u2 - s.*z3 - 73*z2.*z3;
    dz4 = -s.*z4 + 35.2*z1.*z2 - 51.3*z4.*z5;
    dz5 = -s.*z5 + 219*z2.*z3 - 51.3*z4.*z5;
    dz6 = -s.*z6 + 102.6*z4.*z5 - 23*z1.*z6.*u3;
    dz7 = -s.*z7 + 46*z1.*z6.*u3;
    dz8 = 5.8*(s.*z1 - u4) - 3.7*u1 - 4.1*u2 ...
        + s.*(23*z4 + 11*z5 + 28*z6 + 35*z7) - 5*u3.^2 - 0.09;

    DZ = [dz1; dz2; dz3; dz4; dz5; dz6; dz7; dz8];

    ceq = zeros(1, 8*m);
    for i = 1:8
        xdi = w((i-1)*m+1 : i*m);
        ceq((i-1)*m+1 : i*m) = xdi - tf*DZ(i,:);
    end

    c = [];

end
