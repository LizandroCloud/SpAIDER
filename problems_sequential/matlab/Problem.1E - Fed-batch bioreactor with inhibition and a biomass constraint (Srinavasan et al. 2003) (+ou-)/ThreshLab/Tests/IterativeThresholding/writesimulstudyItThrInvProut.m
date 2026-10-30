
if exist('filenameout') ~= 1, filenameout = 'simulstudyItThrInvProut.m'; end
fid = fopen(filenameout,'w');

fprintf(fid,'\nloopnumber = %d;\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n',...
        simul);

fprintf(fid,'vecmseinput = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%15.11f ',vecmseinput(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecmseitdepHTminmse = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%15.11f ',vecmseitdepHTminmse(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecmseitdepHTmingcv = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%15.11f ',vecmseitdepHTmingcv(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecmsefixedHT = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%15.11f ',vecmsefixedHT(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecmseitdepSTminmse = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%15.11f ',vecmseitdepSTminmse(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecmseitdepSTmingcv = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%15.11f ',vecmseitdepSTmingcv(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecmsefixedST = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%15.11f ',vecmsefixedST(simul));
end, fprintf(fid,'];\n');






fprintf(fid,'vecfalsenegitdepHTminmse = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsenegitdepHTminmse(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalsenegitdepHTmingcv = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsenegitdepHTmingcv(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalsenegfixedHT = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsenegfixedHT(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalsenegitdepSTminmse = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsenegitdepSTminmse(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalsenegitdepSTmingcv = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsenegitdepSTmingcv(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalsenegfixedST = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsenegfixedST(simul));
end, fprintf(fid,'];\n');





fprintf(fid,'vecfalsepositdepHTminmse = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsepositdepHTminmse(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalsepositdepHTmingcv = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsepositdepHTmingcv(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalseposfixedHT = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalseposfixedHT(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalsepositdepSTminmse = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsepositdepSTminmse(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalsepositdepSTmingcv = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalsepositdepSTmingcv(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vecfalseposfixedST = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vecfalseposfixedST(simul));
end, fprintf(fid,'];\n');






fprintf(fid,'vectruepositdepHTminmse = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vectruepositdepHTminmse(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vectruepositdepHTmingcv = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vectruepositdepHTmingcv(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vectrueposfixedHT = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vectrueposfixedHT(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vectruepositdepSTminmse = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vectruepositdepSTminmse(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vectruepositdepSTmingcv = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vectruepositdepSTmingcv(simul));
end, fprintf(fid,'];\n');

fprintf(fid,'vectrueposfixedST = [');
for simul=1:nsimul0, fprintf(fid,'\n');
   fprintf(fid,'%d ',vectrueposfixedST(simul));
end, fprintf(fid,'];\n');

fclose(fid);
