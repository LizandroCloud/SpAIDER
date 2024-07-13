function [t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,Wav,adap,frozen,...
    solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,tunn,fref,typ,wfact,reg,est,folder] = spaider_head(optSPAIDER)

    % problem parameters...
    
    t0=optSPAIDER.input.t0;
    tf=optSPAIDER.input.tf;
    is0=optSPAIDER.input.is0;
    isF=optSPAIDER.input.isF;
    nU=optSPAIDER.input.nU;
    x0=optSPAIDER.input.x0;
    u0=optSPAIDER.input.u0;
    u_lb=optSPAIDER.input.u_lb;
    u_ub=optSPAIDER.input.u_ub ;
    
    wname=optSPAIDER.wavelets.wname;
    tptr=optSPAIDER.wavelets.tptr;
    sorh=optSPAIDER.wavelets.sorh ;
    scal=optSPAIDER.wavelets.scal;
    Wav=optSPAIDER.wavelets.procedure ;
    adap=optSPAIDER.wavelets.thresholding;
    
    frozen=optSPAIDER.solver.speed_factor;
    solver=optSPAIDER.solver.name;
    
    plot=optSPAIDER.plot.control;
    plot_state= optSPAIDER.plot.state;
    tf_p=optSPAIDER.plot.tf_p;
    f_tol=optSPAIDER.conv.f_tol;
    f_tol2=optSPAIDER.conv.f_tol2;
    t_ref=optSPAIDER.conv.t_ref;
    tunn=optSPAIDER.wavelets.max_refinement;
    typ=optSPAIDER.problem.tfree;
    wfact=optSPAIDER.wavelets.refinement_position ;
    reg=optSPAIDER.solver.reg;
    est=optSPAIDER.solver.est;
    fref = optSPAIDER.conv.fref;
    folder  = optSPAIDER.folder;
end