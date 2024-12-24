function [s] = write_obj(J, dJ,param)

% write obf1
fid = fopen('obf1.m','w');
fid = fopen('obf1.m','r+');
fprintf(fid,'%10000.s\n', '');
fclose(fid);
fid = fopen('obf1.m','w');
fid = fopen('obf1.m','r+');
fprintf(fid,'%c', 'function [J] = obf1(f,t_ult,t_out,y,stages)');
fprintf(fid,'%c\n', '');
if isempty(param)
else
    for i=1:length(param(:,1))
        fprintf(fid,'%c', param(i,:),';');
        fprintf(fid,'%c\n', '');
    end
end
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'J=',J,';');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'end');

% write obf2
fid = fopen('obf2.m','w');
fid = fopen('obf2.m','r+');
fprintf(fid,'%10000.s\n', '');
fclose(fid);
fid = fopen('obf2.m','w');
fid = fopen('obf2.m','r+');
fprintf(fid,'%c', 'function [J dJ] = obf2(f,t_ult,t_out,df,ws,y,stages)');
fprintf(fid,'%c\n', '');
if isempty(param)
else
    for i=1:length(param(:,1))
        fprintf(fid,'%c', param(i,:),';');
        fprintf(fid,'%c\n', '');
    end
end
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'J=',J,';');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'dJ=',dJ,';');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'end');
s=1;
end