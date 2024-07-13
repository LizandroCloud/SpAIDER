function [SPAIDER] = spaider_main(optSPAIDER,optNLP,optODE)

if optSPAIDER.conv.utarg == 1;
   optSPAIDER.prob.uref=wopt;
   optSPAIDER.prob.tref=tx(:,:,end);
   optSPAIDER.prob.nref=nsx(1,1,end);
else
   optSPAIDER.prob.uref=[1];
    optSPAIDER.prob.tref=[1];
    optSPAIDER.prob.nref=[1];
end

if optSPAIDER.plot.plot_after == 0
if optSPAIDER.casestudy.state==1
    kl=str2double(optSPAIDER.casestudy.iter);
else
    kl=1;
end
for kg = 1:kl
    if optSPAIDER.casestudy.state==1
        disp(sprintf('CASE STUDY: %d', kg));  % writing on the terminal   
    end
if optSPAIDER.wavelets.adap==0
    tptr1= strcat('fix_heuristic_',num2str(optSPAIDER.wavelets.adap),'_2',num2str(optSPAIDER.wavelets.thr));
elseif optSPAIDER.wavelets.adap==1
    tptr1= optSPAIDER.wavelets.tptr;
elseif optSPAIDER.wavelets.adap==2
     tptr1= strcat('norm_based_',num2str(optSPAIDER.wavelets.adap),'_2',num2str(optSPAIDER.wavelets.normfrac));
end

if optSPAIDER.wavelets.Wav==2
     tptr2=strcat('fixed_',int2str(optSPAIDER.prob.is0),'_',int2str(optSPAIDER.prob.isF));
     tptr1= tptr2;
end

file=strcat(tptr1 , '.txt1');
folder = strcat(pwd , '\','results\step_integration_',int2str(optSPAIDER.solver.est),'\regularization_',int2str(optSPAIDER.solver.reg),'\',tptr1,'\',file);
Ax=exist(folder,'file');
if Ax~=0
delete(folder);
end
diary(file);
diary on;
% 
%--------------------------------------------------------------------------
[SPAIDER.wopt, SPAIDER.erropt, SPAIDER.grad,SPAIDER.hessian,SPAIDER.wx,SPAIDER.tx,SPAIDER.nsx,SPAIDER.time,...
    SPAIDER.time_ac,SPAIDER.Fobj,SPAIDER.thrs,SPAIDER.conv,SPAIDER.convv,SPAIDER.detais0,SPAIDER.detais1,...
    SPAIDER.uth,SPAIDER.ld0,SPAIDER.ld1,SPAIDER.md0,SPAIDER.md1,SPAIDER.vd0,SPAIDER.vd1,SPAIDER.difU,...
    SPAIDER.difd,SPAIDER.mseU,SPAIDER.iter] = spaider2_(optSPAIDER,optNLP,optODE);
% [SPAIDER] = spaider2_(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,Wav,adap,frozen,solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,optNLP,optODE,tunn,fref,typ,wfact,reg,est,optSPAIDER);
% 
 [fk,drel] = spaider_rel(optSPAIDER,SPAIDER);
 [fk,trel] = threshold_rel(optSPAIDER,SPAIDER);
 diary off;

file1=strcat(pwd , '\', file);
folder = strcat(pwd , '\','results\step_integration_',int2str(optSPAIDER.solver.est),'\regularization_',int2str(optSPAIDER.solver.reg),'\',tptr1,'\');
 [SUCCESS,MESSAGE,MESSAGEID] = mkdir(folder);
[SUCCESS,MESSAGE,MESSAGEID] = movefile (file1,folder);
if SUCCESS == 0
    mkdir(folder);
    [SUCCESS,MESSAGE,MESSAGEID] = movefile (file1,folder);
end

file2=strcat(folder,tptr1,'.mat');
save (file2) ;
file22=strcat(folder,tptr1,'.xls'); 
drel{1,4}=cell2mat({''}); drel{1,5}=cell2mat({tptr1});
trel{1,4}=cell2mat({''}); trel{1,5}=cell2mat({tptr1});
srel1 = xlswrite(file22, drel, 'Results', 'B1'); % printing relatory
srel2 = xlswrite(file22, trel, 'Thresholds', 'B1');

for k=3:SPAIDER.iter-1
[fk,irel] = iter_rel(optSPAIDER,SPAIDER,k-2);
aba = cell2mat({'Iteration',int2str(k-1)});
srel3 = xlswrite(file22, irel, aba, 'B1');
end


SPAIDER.convhist(kg,1:length(SPAIDER.conv)) = SPAIDER.conv;
if optSPAIDER.wavelets.adap == 0
    SPAIDER.thf(kg) = optSPAIDER.wavelets.thr;
    SPAIDER.timenorm(kg) = SPAIDER.time_ac(length((SPAIDER.tx(end,end,:))))/(optSPAIDER.conv.t_ref);
    optSPAIDER.wavelets.thr = SPAIDER.thf(kg)  + optSPAIDER.wavelets.delt;
elseif optSPAIDER.wavelets.adap == 2 
    SPAIDER.thf(kg) = optSPAIDER.wavelets.normfrac;
    SPAIDER.timenorm(kg) = SPAIDER.time_ac(length((SPAIDER.tx(end,end,:))))/(optSPAIDER.conv.t_ref);
    optSPAIDER.wavelets.normfrac = SPAIDER.thf(kg)  + optSPAIDER.wavelets.delt;
end



end
else 
% plot=1;
for i=2:length((nsx(1,1,:)))
    for k=1:nU
        ncz=(nsx(k,:,i));
        nz = ncz(ncz~=0 ); 
        if isempty(nz)
            nsz=1;
        else
            nsz=nz;
        end
        n_est(k,i-1) = nsx(k,1,i);
         wx(:,:,i)=wx(:,:,i);
    end
end

[ns] = plotp(t0,tf,is0,i+1,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,Wav,adap,frozen,solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,wx,tx,n_est,Fobj,time,time_ac,wfact,tunn);

end


delete *.asv;

end