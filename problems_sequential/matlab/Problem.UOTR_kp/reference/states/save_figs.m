% Define o diretório onde os arquivos serão salvos
saveDir = 'C:\Dropbox\Spaider_v4 _uotr\problems_sequential\matlab\Problem.UOTR_kp\states\state_wv_1';  % Altere conforme necessário
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% Obtém os identificadores (handles) de todas as figuras abertas
figHandles = findall(0, 'Type', 'figure');

% Loop para salvar cada figura
for i = 1:length(figHandles)
    fig = figHandles(i);
    
    % Define nome base dos arquivos
    figName = fullfile(saveDir, sprintf('Figura_%d', i));
    
    % Salva como .fig (formato nativo do MATLAB)
    savefig(fig, [figName '.fig']);
    
    % Salva como .png (ou outro formato)
    exportgraphics(fig, [figName '.png'], 'Resolution', 300);  % Alta qualidade
    
    % Alternativa mais antiga (para versões antes do R2020a):
    % print(fig, [figName '.png'], '-dpng', '-r300');
end

disp('Todas as figuras foram salvas.');
