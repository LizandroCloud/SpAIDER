function [SPAIDER] = spaider_main(in_SPAIDER,optNLP,optODE,SPAIDER)

[ dx y par e_var a_var c_var J cr_eq cr_deq in_SPAIDER.input] = input_model; % routine where are written the model and optimization

in_SPAIDER.solver.est=0;  in_SPAIDER.solver.reg=0;

if strcmp(in_SPAIDER.problem.interpolation,'stage')
    in_SPAIDER.solver.est=1;
else
    in_SPAIDER.solver.est=0;
end

if strcmp(in_SPAIDER.problem.interpolation,'reg_const')
    in_SPAIDER.solver.reg=1;
else
    in_SPAIDER.solver.reg=0;
end

if strcmp(in_SPAIDER.plot.control,'on')
    in_SPAIDER.plot.control=1;
else
    in_SPAIDER.plot.control=0;
end
    
if in_SPAIDER.conv.utarg == 1;
%    in_SPAIDER.prob.uref=wopt;
%    in_SPAIDER.prob.tref=tx(:,:,end);
%    in_SPAIDER.prob.nref=nsx(1,1,end);
else
   in_SPAIDER.prob.uref=[1 1];
    in_SPAIDER.prob.tref=[0 1];
    in_SPAIDER.prob.nref=[2];
end

if strcmp(in_SPAIDER.plot.plot_after,'off')
if strcmp(in_SPAIDER.casestudy.status,'on')
    kl=str2double(in_SPAIDER.casestudy.iter);
else
    kl=1;
end
for kg = 1:kl
    if strcmp(in_SPAIDER.casestudy.status,'on')
        disp(sprintf('CASE STUDY: %d', kg));  % writing on the terminal   
    end

    if strcmp(in_SPAIDER.wavelets.thresholding,'fixed')
        tptr1= strcat('fix_heuristic_',num2str(in_SPAIDER.wavelets.fix_threshold));
    elseif strcmp(in_SPAIDER.wavelets.thresholding,'automatic')
        tptr1= strcat(in_SPAIDER.wavelets.wname,'_',in_SPAIDER.wavelets.sorh,'_',in_SPAIDER.wavelets.scal,'_',in_SPAIDER.wavelets.tptr);
    elseif strcmp(in_SPAIDER.wavelets.thresholding,'norm-based')
        tptr1= strcat('norm_based_',num2str(in_SPAIDER.wavelets.norm_threshold));
    end

    if strcmp(in_SPAIDER.wavelets.procedure,'off')
        tptr2=strcat('fixed_',int2str(in_SPAIDER.input.is0),'_',int2str(in_SPAIDER.input.isF));
        tptr1= tptr2;
    end

    file=strcat(tptr1 , '.txt1');
    folder = strcat(pwd , '\','results\_',in_SPAIDER.problem.interpolation,'\',tptr1,'\');
    Ax=exist(folder,'file');
    
    if Ax~=0
        try
        rmdir(folder);
        catch ME
     end
    end

in_SPAIDER.folder = folder;
[SUCCESS,MESSAGE,MESSAGEID] = mkdir(folder);
if SUCCESS == 0
    mkdir(folder1);
end
diary(file);
diary on;
% 



%--------------------------------------------------------------------------
[SPAIDER.wopt, SPAIDER.ts,SPAIDER.erropt, SPAIDER.grad,SPAIDER.hessian,SPAIDER.wx,SPAIDER.tx,SPAIDER.nsx,SPAIDER.time,...
    SPAIDER.time_ac,SPAIDER.Fobj,SPAIDER.thrs,SPAIDER.conv,SPAIDER.convv,SPAIDER.detais0,SPAIDER.detais1,...
    SPAIDER.uth,SPAIDER.ld0,SPAIDER.ld1,SPAIDER.md0,SPAIDER.md1,SPAIDER.vd0,SPAIDER.vd1,SPAIDER.difU,...
    SPAIDER.difd,SPAIDER.mseU,SPAIDER.iter,SPAIDER.condh,SPAIDER.time_state,SPAIDER.state,SPAIDER.psnr] = spaider(in_SPAIDER,optNLP,optODE);
% [SPAIDER] = spaider2_(t0,tf,is0,isF,nU,x0,u0,u_lb,u_ub,wname,tptr,sorh,scal,Wav,adap,frozen,solver,plot,plot_state,tf_p,f_tol,f_tol2,t_ref,optNLP,optODE,tunn,fref,typ,wfact,reg,est,in_SPAIDER);
% 
 if in_SPAIDER.problem.tfree == 1
     in_SPAIDER.input.nU=in_SPAIDER.input.nU-0;
 end

 if strcmp(in_SPAIDER.problem.type,'parameter estimation')

 else

 [fk,drel] = spaider_rel(in_SPAIDER,SPAIDER);
  if strcmp(in_SPAIDER.wavelets.procedure,'on')
      
      if strcmp(in_SPAIDER.wavelets.rel,'on')
 [fk,trel] = threshold_rel(in_SPAIDER,SPAIDER);
      end
  end
 diary off;

% file1=strcat(pwd , '\', file);
% folder = strcat(pwd , '\','results\step_integration_',int2str(in_SPAIDER.solver.est),int2str(in_SPAIDER.solver.reg),'\',tptr1,'\');
 [SUCCESS,MESSAGE,MESSAGEID] = mkdir(folder);
[SUCCESS,MESSAGE,MESSAGEID] = movefile (file,folder);
if SUCCESS == 0
    mkdir(folder);
    [SUCCESS,MESSAGE,MESSAGEID] = movefile (file,folder);
end

%  movefile ('*.xls' , folder);
file2=strcat(folder,tptr1,'.mat');
save (file2) ;
file22=strcat(folder,tptr1,'.xls'); 
file33 = strcat(tptr1,'.xls'); 
drel{1,4}=cell2mat({''}); drel{1,5}=cell2mat({tptr1});
trel{1,4}=cell2mat({''}); trel{1,5}=cell2mat({tptr1});
srel1 = xlswrite(file33, drel, 'Results', 'B1'); % printing relatory
% movefile (file33 , folder)
 if strcmp(in_SPAIDER.wavelets.rel,'on')
srel2 = xlswrite(file33, trel, 'Thresholds', 'B1');
 end
 
for k=3:SPAIDER.iter-1
[fk,irel] = iter_rel(in_SPAIDER,SPAIDER,k-2);
aba = cell2mat({'Iteration',int2str(k-1)});
srel3 = xlswrite(file33, irel, aba, 'B1');
end

 if strcmp(in_SPAIDER.wavelets.procedure,'on') 
     try
    movefile (file33 , folder)
    
     catch ME
     end
 end

 if strcmp(in_SPAIDER.wavelets.rel,'on')
     try
    movefile ('*.xls' , folder)
    
     catch ME
     end
 end

end

SPAIDER.convhist(kg,1:length(SPAIDER.conv)) = SPAIDER.conv;
if strcmp(in_SPAIDER.wavelets.thresholding,'fixed')
    SPAIDER.thf(kg) = in_SPAIDER.wavelets.fix_threshold;
    SPAIDER.timenorm(kg) = SPAIDER.time_ac(length((SPAIDER.tx(end,end,:))))/(in_SPAIDER.conv.t_ref);
    in_SPAIDER.wavelets.thr = SPAIDER.thf(kg)  + in_SPAIDER.wavelets.delt;
elseif strcmp(in_SPAIDER.wavelets.thresholding,'norm-based')
    SPAIDER.thf(kg) = in_SPAIDER.wavelets.norm_threshold;
    SPAIDER.timenorm(kg) = SPAIDER.time_ac(length((SPAIDER.tx(end,end,:))))/(in_SPAIDER.conv.t_ref);
    in_SPAIDER.wavelets.normfrac = SPAIDER.thf(kg)  + in_SPAIDER.wavelets.delt;
end

    if in_SPAIDER.problem.tfree == 1
     in_SPAIDER.input.nU=in_SPAIDER.input.nU+1;
    end

end
else 
%  load SPAIDER;
in_SPAIDER.plot.control = 1;
% plot=1;
for i=2:length((SPAIDER.nsx(1,1,:)))
    for k=1:in_SPAIDER.input.nU
        ncz=(SPAIDER.nsx(k,:,i));
        nz = ncz(ncz~=0 ); 
        if isempty(nz)
            nsz=1;
        else
            nsz=nz;
        end
        n_est(k,i-1) = SPAIDER.nsx(k,1,i);
        wxd(:,:,i)=SPAIDER.wx(:,:,i);
    end
    

end
% [ns]     = plotp(in_SPAIDER,optNLP,optODE);
% in_SPAIDER.plot.tf_p =0.2;
[ns] = plotp(in_SPAIDER.input.t0,in_SPAIDER.input.tf,in_SPAIDER.input.is0,i+1,in_SPAIDER.input.nU,in_SPAIDER.input.x0,...
in_SPAIDER.input.u0,in_SPAIDER.input.u_lb,in_SPAIDER.input.u_ub,...
in_SPAIDER.wavelets.wname, in_SPAIDER.wavelets.tptr, in_SPAIDER.wavelets.sorh ,...
in_SPAIDER.wavelets.scal ,in_SPAIDER.wavelets.procedure, in_SPAIDER.wavelets.thresholding,...
in_SPAIDER.solver.speed_factor, in_SPAIDER.solver.name,in_SPAIDER.plot.control,in_SPAIDER.plot.state,...
in_SPAIDER.plot.tf_p, in_SPAIDER.conv.f_tol, in_SPAIDER.conv.f_tol2,...
in_SPAIDER.conv.t_ref,wxd,SPAIDER.tx,n_est,SPAIDER.Fobj,SPAIDER.time,  SPAIDER.time_ac,...
in_SPAIDER.wavelets.refinement_position,  in_SPAIDER.wavelets.max_refinement,in_SPAIDER);

end


delete *.asv;

end