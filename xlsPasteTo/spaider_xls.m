% Create some arbitrary graphics
% f1 = figure; peaks; f2 = figure; membrane;
f1 = figure(1); peaks; set (1,'Visible','off');

% Connect to Excel, make it visible and add a worksheet
xl = actxserver('Excel.Application'); set(xl,'Visible',0);
xl.Workbooks.Add(1); 
xls = xl.ActiveSheet;
is='1';
% Paste in the MATLAB figures
print(1, '-dbitmap');  set (1,'Visible','off');
xls.Range('E3').PasteSpecial;
xls.SaveAs([pwd '\iteration',is,'.xls']);
invoke(xls, 'Saveas', [pwd '\iteration',is,'.xls']);
invoke(Excel, 'Quit');
% print(f2, '-dbitmap'); xls.Range('I3').PasteSpecial;
% xl.SaveAs(pwd,'example2');
% xls.Shapes.Item(1).PictureFormat.CropLeft  = 30;
% xls.Shapes.Item(1).PictureFormat.CropRight  = 30;
% xls.Shapes.Item(1).Height  = 200;
% xls.Shapes.Item(1).Left = xls.Range('E3').Left;
