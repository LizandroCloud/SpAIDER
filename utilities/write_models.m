function [s] = write_models(in_SPAIDER,eq, eq_sens, param, n_p, n_s, n_c, n_a, y)
% write modelo_eq (for numerical integration)
fid = fopen('modelo_eq.m','w');
fid = fopen('modelo_eq.m','r+');
fprintf(fid,'%10000.s\n', '');
fclose(fid);
fid = fopen('modelo_eq.m','r+');
fprintf(fid,'%c', 'function [ dx ] = modelo_eq( t_out, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,tal,est,yi,utarg,tstarg,nstarg)');
fprintf(fid,'%c\n', '');
switch (in_SPAIDER.problem.interpolation) 
     case 'reg_const'
        fprintf(fid,'%c', '[u] = ver_reg_const(ws,t,te,ks,ns,reg,est);');
     case 'linear'
        fprintf(fid,'%c', '[u] = ver_reg_est0(ws,yi,t_out);');    
     case 'spline'   
        fprintf(fid,'%c', '[u] = ver_reg_est0(ws,yi,t_out);');    
     case 'nearest'   
        fprintf(fid,'%c', '[u] = ver_reg_est0(ws,yi,t_out);');    
     case 'stage'   
        fprintf(fid,'%c', '[u] = ver_reg_est1(ws,yi,t_out);');    
end         
fprintf(fid,'%c\n', '');
switch (in_SPAIDER.problem.interpolation) 
     case 'reg_const'
        fprintf(fid,'%c', '[uref] = ver_reg_const(ws,t,te,ks,ns,reg,est);');
     case 'linear'
        fprintf(fid,'%c', '[uref] = ver_reg_est0(ws,yi,t_out);');    
     case 'spline'   
        fprintf(fid,'%c', '[uref] = ver_reg_est0(ws,yi,t_out);');    
     case 'nearest'   
        fprintf(fid,'%c', '[uref] = ver_reg_est0(ws,yi,t_out);');    
     case 'stage'   
        fprintf(fid,'%c', '[uref] = ver_reg_est1(ws,yi,t_out);');    
end 
fprintf(fid,'%c\n', '');
for i=1:n_p
    fprintf(fid,'%c', param(i,:),';');
    fprintf(fid,'%c\n', '');
end
fprintf(fid,'%c\n', '');
if n_a==0
else
    for i=1:n_a
         a=int2str(i);
         fprintf(fid,'%c', 'y(',a,')=','(',y(i,:),')',';');
         fprintf(fid,'%c\n', '');
    end
end
fprintf(fid,'%c\n', '');
for i=1:n_s
    a=int2str(i);
    fprintf(fid,'%c', 'dx(',a,')=','(',eq(i,:),')*tal',';');
    fprintf(fid,'%c\n', '');
end
fprintf(fid,'%c','dx=transp(dx);');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'end');
fclose(fid);
s=1;
% write modelo_eq (for numerical integration and sensitivity)
fid = fopen('modelo_eq_sens.m','w');
fid = fopen('modelo_eq_sens.m','r+');
fprintf(fid,'%10000.s\n', '');
fclose(fid);
fid = fopen('modelo_eq_sens.m','w');
fid = fopen('modelo_eq_sens.m','r+');
fprintf(fid,'%c', 'function [ dx ] = modelo_eq_sens( t_out, x, ws, ks,ts,reg,ns,ne,te,n_s,n_c,ns_a,ts_frozen,points,uws,tal,est,yi,utarg,tstarg,nstarg)');
fprintf(fid,'%c\n', '');
switch (in_SPAIDER.problem.interpolation) 
     case 'reg_const'
        fprintf(fid,'%c', '[u] = ver_reg_const(ws,t,te,ks,ns,reg,est);');
     case 'linear'
        fprintf(fid,'%c', '[u] = ver_reg_est0(ws,yi,t_out);');    
     case 'spline'   
        fprintf(fid,'%c', '[u] = ver_reg_est0(ws,yi,t_out);');    
     case 'nearest'   
        fprintf(fid,'%c', '[u] = ver_reg_est0(ws,yi,t_out);');    
     case 'stage'   
        fprintf(fid,'%c', '[u] = ver_reg_est1(ws,yi,t_out);');    
end 
fprintf(fid,'%c\n', '');
switch (in_SPAIDER.problem.interpolation) 
     case 'reg_const'
        fprintf(fid,'%c', '[uref] = ver_reg_const(ws,t,te,ks,ns,reg,est);');
     case 'linear'
        fprintf(fid,'%c', '[uref] = ver_reg_est0(ws,yi,t_out);');    
     case 'spline'   
        fprintf(fid,'%c', '[uref] = ver_reg_est0(ws,yi,t_out);');    
     case 'nearest'   
        fprintf(fid,'%c', '[uref] = ver_reg_est0(ws,yi,t_out);');    
     case 'stage'   
        fprintf(fid,'%c', '[uref] = ver_reg_est1(ws,yi,t_out);');    
end 
fprintf(fid,'%c\n', '');
for i=1:n_p
    fprintf(fid,'%c', param(i,:),';');
    fprintf(fid,'%c\n', '');
end
fprintf(fid,'%c\n', '');
if n_a==0
else
for i=1:n_a
    a=int2str(i);
    fprintf(fid,'%c', 'y(',a,')=','(',y(i,:),')',';');
    fprintf(fid,'%c\n', '');
end
end
fprintf(fid,'%c\n', '');
for i=1:n_s
    a=int2str(i);
    fprintf(fid,'%c', 'dx(',a,')=','(',eq(i,:),')*tal',';');
    fprintf(fid,'%c\n', '');
end
fprintf(fid,'%c','dx=transp(dx);');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','jc1 = [');
fprintf(fid,'%c\n', '');
for i=1:n_s
    end1=(length(eq_sens(i,:)-3));
    fprintf(fid,'%c',eq_sens(i,:));
    fprintf(fid,'%c\n', '');
end
fprintf(fid,'%c',']*tal;');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','jk=double(jc1);');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','j0=0;');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','jk(:,1:n_s)=jc1(:,1:n_s);');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','jk(:,n_s+1:n_s+n_c)=0;'); 
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','for j=1:n_c');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','for is = 1:ns_a(j)');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','jk(:,n_s+j)=uws(j,is)*jc1(:,n_s+j);');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','k=1:(n_s);');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','dx(n_s*is+k+j0) = jk(k,1:n_s)*x(n_s*is+1+j0:n_s*is+n_s+j0)+ jk(k,n_s+j);');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','end ');                
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','j0=length(dx)-n_s;');      
fprintf(fid,'%c\n', '');
fprintf(fid,'%c','end '); 
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'end');
fclose(fid);
s=1;
end