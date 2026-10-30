
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulation study Iterative Thresholding for Inverse Problems %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Investigates and compares the following methods:
%  1) it.dep. ST with argmin(GCVsoft) in each iteration
%  2) it.dep. HT with argmin(GCVhard) in each iteration
%  3) it.dep. ST with argmin(MSEsoft) in each iteration
%  4) it.dep. HT with argmin(MSEhard) in each iteration
%  5) fixed ST with min MSE threshold
%  6) fixed HT with min MSE threshold
%
% typical behaviors:
%  1) large, irregular oscillatory MSE
%  2) smooth, almost monotone decay of MSE
%  3) small, regular oscillatory MSE
%  4) mostly monotone decay of MSE with sudden jumps (change of attraction
%     point in non-convex optimization problem?)
%
% Each method is reported on 
%  a) limiting MSE
%  b) limiting number of correct positives (true discoveries)
%  c) limiting number of false negatives (missing discoveries)
%  d) limiting number of false positives (false discoveries)


nsimul = 200;
realrandom = 1;
% nsimul0 is the number of simulation loops that is read from filenameout (if
% it exists already) and it is also the number of loops that will be written.
% nsimul0 may be larger than nsimul (but not smaller) if the existing file
% contains more simulations than the new simulation study will cover. This way,
% old simulations are read and re-written without being lost, simply because
% the current run has less loops.

filenameout = 'simulstudyItThrInvProut.m';

if exist(filenameout) ~= 2,
   nsimul0 = 0;
   initsimulstudyItThrInvPr
   nsimul0 = nsimul;
   simul = 0;
   writesimulstudyItThrInvProut
else
   simulstudyItThrInvProut
   nsimul0 = length(vecmseinput);
   if nsimul0 < nsimul,
      initsimulstudyItThrInvPr
      nsimul0 = nsimul;
   end
end

for simul = 174:nsimul,
   disp('---------------------------------')
   disp(['start simulation ' num2str(simul)])
   disp('---------------------------------')

   randn('state',10*simul);
   rand('state',10*simul);

   simulationItThrInvPr

   summarizesimulItThrInvPr
   writesimulstudyItThrInvProut
end
