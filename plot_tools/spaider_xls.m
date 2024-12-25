function [pic2xls] = spaider_xls(ixs,is,nU,folder)

% Create some arbitrary graphics
% f1 = figure; peaks; f2 = figure; membrane;
%  figure(1); peaks; 
%  figure(2); membrane; 
%  set (1,'Visible','off');  
%  set (2,'Visible','off');

% Connect to Excel, make it visible and add a worksheet
ME.message = ''; 
obs=0;
xl = actxserver('Excel.Application'); set(xl,'Visible',1);
xl.Workbooks.Add(1); 
xls = xl.ActiveSheet;
% is='5';
% Paste in the MATLAB figures
j = 1;
for i=1:length(ixs)
print(ixs(i), '-dbitmap'); 
position = strcat('B',int2str(j));
xls.Range(position).PasteSpecial;
j = j+30;
end
set (ixs,'Visible','off');
iss = num2str(is-3);
nU1 = num2str(nU);
filename = strcat('iteration',iss,'control',nU1,'.xlsx');
xls.SaveAs([pwd '\' filename]);
invoke(xl, 'Quit');
% invoke(xls, 'Saveas', [pwd '\iteration',is,'.xls']);
%   while obs == 0
%  
%   try
%       
%     movefile (filename , folder);
%       catch ME
%                
%   end
%    obs= strcmp(ME.message,'');
%    pause
%   end
% movefile (filename , folder)
% print(f2, '-dbitmap'); xls.Range('I3').PasteSpecial;
% xl.SaveAs(pwd,'example2');
% xls.Shapes.Item(1).PictureFormat.CropLeft  = 30;
% xls.Shapes.Item(1).PictureFormat.CropRight  = 30;
% xls.Shapes.Item(1).Height  = 200;
% xls.Shapes.Item(1).Left = xls.Range('E3').Left
savefigs(ixs,is,nU);
set (ixs,'Visible','off');
pic2xls = 1;
% movefile ('*.png' , folder);
% movefile ('*.emf' , folder);
% movefile ('*.eps' , folder);
end
