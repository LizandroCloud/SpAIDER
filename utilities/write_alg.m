function [s] = write_alg(eq, eq_sens, param, n_p, n_s, n_c, n_a, y)

% write obf1  % numerical derivatives
fid = fopen('alg_eq.m','w');
fid = fopen('alg_eq.m','r+');
fprintf(fid,'%10000.s\n', '');
fclose(fid);
fid = fopen('alg_eq.m','r+');
fprintf(fid,'%c', 'function [yout] = alg_eq(x,u)');
fprintf(fid,'%c\n', '');
if isempty(param)
else
for i=1:length(param(:,1))
    fprintf(fid,'%c', param(i,:),';');
    fprintf(fid,'%c\n', '');
end
end

% writting algebraic equations...
fprintf(fid,'%c\n', '');
if n_a==0
         fprintf(fid,'%c\n', '');
         fprintf(fid,'%c', 'y(1)=1;');
         fprintf(fid,'%c\n', '');
else
    for i=1:n_a
         a=int2str(i);
         fprintf(fid,'%c\n', '');
         fprintf(fid,'%c', 'y(',a,')=','(',y(i,:),')',';');
         fprintf(fid,'%c\n', '');
    end
end
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'yout=y;');
fprintf(fid,'%c\n', '');
fprintf(fid,'%c', 'end');
s=1;
end