disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('  Welcome to Spaider(Sparse Algorithms Integrated to Differential Equations Routines), version 1.0.1, 8 November 2012   ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  
disp('             Copyright (c) Lizandro de Sousa Santos')
disp('Currently running the startup file, located in Spaider/Problem_name/startup.m')
disp('...')
P1 = pwd; 
cd .. 
cd ..
cd ..
P2= pwd;
path0 = P2;
path(path,path0);
path(path,[path0 '/dasslc/']);
path(path,[path0 '/odes/']);
path(path,[path0 '/plot_tools/']);
path(path,[path0 '/test_routines']);
path(path,[path0 '/utilities']);
path(path,[path0 '/wavelets']);
path(path,[path0 '/functions']);

disp('---------------------------')
disp('Spaider path was configurated')
disp('---------------------------')
disp('Your current working directory/folder is:')
disp(pwd)
disp(P1)
cd (P1);
