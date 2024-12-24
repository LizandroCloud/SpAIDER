function [tex tex_o tex_1] = f_te(i,te,ks)
tex = (te(i,:));
tex_o = te(i,ks);
tex_1 = te(i,ks+1);
end