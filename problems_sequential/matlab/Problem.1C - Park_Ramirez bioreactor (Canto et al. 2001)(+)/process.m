function sys = process(t,x,flag,u,p)
switch flag
    case 0
        % Equações dinâmicas
        x1 = x(1);
        x2 = x(2);
        x3 = x(3);
        x4 = x(4);
        x5 = x(5);
        u1 = u(1);

        % Cálculo das derivadas
        dx = zeros(5,1);
        dx(1) = (4.75*(21.88*x4/( (x4 + 0.4)*(x4+62.5) ) ) / (0.12 + (21.88*x4/( (x4 + 0.4)*(x4+62.5) ) )))*(x2-x1)-u1*x1/x5;
        dx(2) = ( ( x4/(0.1+x4) )*exp(-5.01*x4) )*x3 - u1*x2/x5;
        dx(3) = (21.88*x4/( (x4 + 0.4)*(x4+62.5) ) )*x3 - u1*x3/x5;
        dx(4) = -7.3*(21.88*x4/( (x4 + 0.4)*(x4+62.5) ) )*x3 + u1*(20-x4)/x5;
        dx(5) = u1;
        sys = dx;

    case 5
        sys = [0  0  1  5  1];
    otherwise
        sys = [];
end
