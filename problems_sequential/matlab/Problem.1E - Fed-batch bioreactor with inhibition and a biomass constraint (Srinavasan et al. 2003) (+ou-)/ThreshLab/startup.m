disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
disp('  Welcome to ThreshLab, version 3.2.1, 1 March 2012   ')
disp('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')  
disp('                    Copyright (c) Maarten Jansen')
disp('Currently running the ThreshLab startup file, located in ThreshLab/startup.m')
disp('...')
% cd <yourmatlabpath>/ThreshLab/
% uncomment the line above if you want to run startup from another directory
% (make sure that <yourmatlabpath>/ThreshLab/ is already in your path)
path0 = pwd;
path(path,path0);
path(path,[path0 '/General']);
path(path,[path0 '/Plots']);
path(path,[path0 '/WT/']);
path(path,[path0 '/Noise']);
path(path,[path0 '/Threshold']);
path(path,[path0 '/Threshold/GCV']);
path(path,[path0 '/Threshold/MSE']);
path(path,[path0 '/Threshold/Risk']);
path(path,[path0 '/Threshold/FDR']);
path(path,[path0 '/Threshold/EBayesThresh']);
path(path,[path0 '/Threshold/IterativeThresholding']);
path(path,[path0 '/MRF']);
path(path,[path0 '/Poisson']);
path(path,[path0 '/Poisson/PoissonTests']);
path(path,[path0 '/Tests']);
path(path,[path0 '/Tests/GCVdemo']);
path(path,[path0 '/Tests/Grids']);
path(path,[path0 '/Tests/minMSEforHT']);
path(path,[path0 '/Tests/IterativeThresholding']);
path(path,[path0 '/OldNames']);
path(path,[path0 '/ThreshLab2']);
path(path,[path0 '/ThreshLab2/GeneralGeometry']);
path(path,[path0 '/ThreshLab2/Triangulations']);
path(path,[path0 '/ThreshLab2/Polylines']);
path(path,[path0 '/ThreshLab2/Irregular']);
path(path,[path0 '/ThreshLab2/Irregular/1D']);
path(path,[path0 '/ThreshLab2/Irregular/2D']);
path(path,[path0 '/ThreshLab2/MultiscaleKernel']);
path(path,[path0 '/ThreshLab2/NormalOffsets']);
path(path,[path0 '/ThreshLab2/NormalOffsets/NM1D']);
path(path,[path0 '/ThreshLab2/NormalOffsets/NMpixel']);
path(path,[path0 '/ThreshLab2/NormalOffsets/NMupdate1']);
path(path,[path0 '/ThreshLab2/NormalOffsets']);
path(path,[path0 '/ThreshLab2/ContinuousUnbalancedHaar']);
disp('---------------------------')
disp('ThreshLab path has been set')
disp('---------------------------')
disp('Your current working directory/folder is:')
disp(pwd)
