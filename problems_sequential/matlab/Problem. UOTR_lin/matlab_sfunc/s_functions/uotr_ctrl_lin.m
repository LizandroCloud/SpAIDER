function dxdt=uotr_ctrl_lin(t,x,u,dFdx,dFdu,x0,u0)
    dxdt=dFdx*(x-x0)+dFdu*(u-u0);
end