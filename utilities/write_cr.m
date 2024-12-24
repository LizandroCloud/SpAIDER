function [s] = write_cr(cr_eq, d_cr_eq, cr_deq, d_cr_deq, param, cr_eq_type ,cr_deq_type)

% write obf1  % numerical derivatives
fid = fopen('cr1.m','w');
fid = fopen('cr1.m','r+');
fprintf(fid,'%10000.s\n', '');
fclose(fid);
fid = fopen('cr1.m','r+');
fprintf(fid,'%c', 'function [ceq c] = cr1(f,t_ult,t_out,u,y,stages)');
fprintf(fid,'%c\n', '');
if isempty(param)
else
for i=1:length(param(:,1))
    fprintf(fid,'%c', param(i,:),';');
    fprintf(fid,'%c\n', '');
end
end
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'timeL=0',';');  % initial time counter
fprintf(fid,'%c\n', '');
    % equality
    fprintf(fid,'%c\n', '');
    fprintf(fid,'%c', '%equality');
    fprintf(fid,'%c\n', ''); 
    k=0;
    for i=1:length(cr_eq(:,1)) % for each eq constraint
        a=int2str(i);
        if cr_eq_type(i)==1    % not null
            fprintf(fid,'%c', 'ceq(',a,') = ',cr_eq(i,:),';');
            fprintf(fid,'%c\n', '');
            k=k+1;
        elseif cr_eq_type(i)==99  % null
            fprintf(fid,'%c', 'ceq = ',cr_eq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', '');
            k=k+1;
        else  % looping
            fprintf(fid,'%c', 'for time=1+0+',int2str(k),':length(f(:,1))+0+',int2str(k));
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'ceq(time+timeL) = ',cr_eq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'end');
            fprintf(fid,'%c\n', '')
            if i<length(cr_eq(:,1))
                fprintf(fid,'%c', 'timeL=length(f(:,1));',';');
                fprintf(fid,'%c\n', '');
            else
%                 fprintf(fid,'%c\n', '');
%                 fprintf(fid,'%c', 'ceq = spline([1:length(ceq)],ceq,1:300);');  
                fprintf(fid,'%c\n', '');
            end
        end
    end

    fprintf(fid,'%c', 'timeL=0',';');  % initial time counter
    fprintf(fid,'%c\n', '');
    
    % inequality
    fprintf(fid,'%c\n', '');
    fprintf(fid,'%c', '%inequality');
    fprintf(fid,'%c\n', ''); 
    k=0;
    
    for i=1:length(cr_deq(:,1)) % for each eq constraint
        a=int2str(i);
        if cr_deq_type(i)==1 % not null
            fprintf(fid,'%c', 'c(',a,') = ',cr_deq(i,:),';');
            fprintf(fid,'%c\n', '');
            k=k+1;
        elseif cr_deq_type(i)==99 %  null
            fprintf(fid,'%c', 'c = ',cr_deq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', '');
            k=k+1;
        else % looping
            fprintf(fid,'%c', 'for time=1+0+',int2str(k),':length(f(:,1))+0+',int2str(k));
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'c(time+timeL) = ',cr_deq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'end');
            fprintf(fid,'%c\n', '');
            if i<length(cr_deq(:,1))
                    fprintf(fid,'%c', 'timeL=length(f(:,1));',';');
                    fprintf(fid,'%c\n', '');
            else
%                 fprintf(fid,'%c\n', '');
%                 fprintf(fid,'%c', 'c = spline([1:length(c)],c,1:300);');  
                fprintf(fid,'%c\n', '');
            end
        end
    end
      
    fprintf(fid,'%c\n', '');
    fprintf(fid,'%c\n', '');
    fprintf(fid,'%c', 'end');
    
% write obf2 % analitical derivatives
fid = fopen('cr2.m','w');
fid = fopen('cr2.m','r+');
fprintf(fid,'%10000.s\n', '');
fclose(fid);
fid = fopen('cr2.m','r+');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'function [ceq c dceq dc] = cr2(f,t_ult,t_out,df,u,y,stages)');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c\n', '');
if isempty(param)
else
for i=1:length(param(:,1))
    fprintf(fid,'%c', param(i,:),';');
    fprintf(fid,'%c\n', '');
end
end
fprintf(fid,'%c\n', '');

fprintf(fid,'%c', 'timeL=0',';');  % initial time counter
fprintf(fid,'%c\n', '');
    
    % equality
    fprintf(fid,'%c\n', '');
    fprintf(fid,'%c', '%equality');
    fprintf(fid,'%c\n', ''); 
    k=0;
    for i=1:length(cr_eq(:,1)) % for each eq constraint
        a=int2str(i);
        if cr_eq_type(i)==1
            fprintf(fid,'%c', 'ceq(',a,') = ',cr_eq(i,:),';');
            fprintf(fid,'%c\n', '');
            k=k+1;
        elseif cr_eq_type(i)==99
            fprintf(fid,'%c', 'ceq = ',cr_eq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', '');
            k=k+1;
        else
            fprintf(fid,'%c', 'for time=1+0+',int2str(k),':length(f(:,1))+0+',int2str(k));
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'ceq(time+timeL) = ',cr_eq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'end');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', '');
            if i<length(cr_eq(:,1))
                fprintf(fid,'%c', 'timeL=length(f(:,1));');
                fprintf(fid,'%c\n', '');
                fprintf(fid,'%c\n', '');
            else
%                 fprintf(fid,'%c\n', '');
%                 fprintf(fid,'%c', 'ceq = spline([1:length(ceq)],ceq,1:300);');  
                fprintf(fid,'%c\n', '');
            end
        end
    end
    fprintf(fid,'%c', 'timeL=0',';');  % initial time counter
    fprintf(fid,'%c\n', '');
    
    % inequality
    fprintf(fid,'%c\n', '');
    fprintf(fid,'%c', '%inequality');
    fprintf(fid,'%c\n', '');
    k=0;
    for i=1:length(cr_deq(:,1)) % for each eq constraint
        a=int2str(i);
        if cr_deq_type(i)==1
            fprintf(fid,'%c', 'c(',a,') = ',cr_deq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', '');
            k=k+1;
        elseif cr_deq_type(i)==99
            fprintf(fid,'%c', 'c = ',cr_deq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', '');
            k=k+1;
        else
            fprintf(fid,'%c', 'for time=1+0+',int2str(k),':length(f(:,1))+0+',int2str(k));
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'c(time+timeL) = ',cr_deq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'end');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', ''); 
            if i<length(cr_deq(:,1))
                     fprintf(fid,'%c', 'timeL=length(f(:,1));',';');
                     fprintf(fid,'%c\n', '');
                     fprintf(fid,'%c\n', '');
            else
%                     fprintf(fid,'%c\n', '');
%                 fprintf(fid,'%c', 'c = spline([1:length(c)],c,1:300);');  
                fprintf(fid,'%c\n', '');
            end
        end
        fprintf(fid,'%c\n', '');
    end
    % equality  derivative
     fprintf(fid,'%c', 'timeL=0',';');  % initial time counter
    fprintf(fid,'%c\n', '');
    fprintf(fid,'%c', '%equality  derivative');
    fprintf(fid,'%c\n', '');
    k=0;
    for i=1:length(cr_eq(:,1)) % for each eq constraint
        a=int2str(i);
        if cr_eq_type(i)==1
            fprintf(fid,'%c', 'dceq(',a,',:) = ',d_cr_eq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', '');
            k=k+1;
        elseif cr_eq_type(i)==99
            fprintf(fid,'%c', 'dceq = ',d_cr_eq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', '');
            k=k+1;
        else
            fprintf(fid,'%c', 'for time=1+0+',int2str(k),':length(f(:,1))+0+',int2str(k));
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'dceq(time+timeL,:) = ',d_cr_eq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'end');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c\n', '');
            if i<length(cr_eq(:,1))
                fprintf(fid,'%c', 'timeL=length(f(:,1));');
                fprintf(fid,'%c\n', '');
                fprintf(fid,'%c\n', '');
            else
                fprintf(fid,'%c\n', '');
                fprintf(fid,'%c\n', '');
                fprintf(fid,'%c', 'for i=1:length(df(1,1,:))');
                fprintf(fid,'%c\n', '');
%                 fprintf(fid,'%c', 'dcceq(:,i) = spline([1:length(dceq(:,1))],dceq(:,i),1:300);');
%                 fprintf(fid,'%c\n', '');
                fprintf(fid,'%c', 'end');
                fprintf(fid,'%c\n', '');
                fprintf(fid,'%c', 'dceq=transp(dcceq);');
                fprintf(fid,'%c\n', '');
            end
        end
    
    fprintf(fid,'%c\n', '');
    end
    fprintf(fid,'%c', 'timeL=0',';');  % initial time counter
    % inequality  derivative
    fprintf(fid,'%c\n', '');
    fprintf(fid,'%c', '%inequality  derivative');
    fprintf(fid,'%c\n', '');
    k=0;
    for i=1:length(cr_deq(:,1)) % for each eq constraint
        a=int2str(i);
        if cr_deq_type(i)==1
            fprintf(fid,'%c', 'dc(',a,',:) = ',d_cr_deq(i,:),';');
            fprintf(fid,'%c\n', '');
            k=k+1;
        elseif cr_deq_type(i)==99
            fprintf(fid,'%c', 'dc = ',d_cr_deq(i,:),';');
            fprintf(fid,'%c\n', '');
            k=k+1;
        else
            fprintf(fid,'%c', 'for time=1+0+',int2str(k),':length(f(:,1))+0+',int2str(k));
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'dc(time+timeL,:) = ',d_cr_deq(i,:),';');
            fprintf(fid,'%c\n', '');
            fprintf(fid,'%c', 'end');
            fprintf(fid,'%c\n', '');
            if i<length(cr_deq(:,1))
                fprintf(fid,'%c', 'timeL=length(f(:,1));');
                fprintf(fid,'%c\n', '');
            else
                fprintf(fid,'%c\n', '');
                fprintf(fid,'%c', 'for i=1:length(df(1,1,:))');
%                 fprintf(fid,'%c\n', '');
%                 fprintf(fid,'%c', 'dcc(:,i) = spline([1:length(dc(:,1))],dc(:,i),1:300);');
                fprintf(fid,'%c\n', '');
                fprintf(fid,'%c', 'end');
                fprintf(fid,'%c\n', '');
            end
        end
    end
 fprintf(fid,'%c\n', '');
%  fprintf(fid,'%c', 'dicc(:,i) = spline([1:length(dc(:,1))],dc(:,i),1:300);');
 fprintf(fid,'%c\n', '');
 fprintf(fid,'%c', 'dc=transp(dc);');
 fprintf(fid,'%c\n', '');
 fprintf(fid,'%c', 'dceq=transp(dceq);');
 fprintf(fid,'%c\n', '');
 fprintf(fid,'%c', 'end');
 s=1;
end