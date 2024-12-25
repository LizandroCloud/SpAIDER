for i=1:length((nsx(1,1,:)))-1
    for k=1:nU
        ncz=(nsx(k,:,i));
        nz = ncz(ncz~=0 ); 
        if isempty(nz)
            nsz=1;
        else
            nsz=nz;
        end
        n_est(k,i) = nsx(k,1,i+1);
         wx(:,:,i)=wx(:,:,i);
    end
end

  [ns] = plotp(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,Wav,adap,frozen,solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,wx,tx,n_est,Fobj,time);
