function [to tf dtt]  = f_dt(ti,ks,disc);

to = ti(ks);
tf = ti(ks+1);
dtt = (tf-to)/1;

end