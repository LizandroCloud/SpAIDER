figHandles = findall(0, 'Type', 'figure'); % Encontra todas as figuras abertas
outputFolder = 'C:\Dropbox\Spaider_v4 _uotr\problems_sequential\matlab\Problem.UOTR_kp\reference'; % Defina o diretório onde deseja salvar

for i = 1:length(figHandles)
    fig = figHandles(i);
    filename = fullfile(outputFolder, sprintf('Figura_%d.png', i)); % Nome do arquivo
    saveas(fig, filename); % Salva como PNG (pode mudar o formato)
end

disp('Todas as figuras foram salvas com sucesso!');
