function convert_to_process(dx_strings, x0)
    % dx_strings: cell array of strings com as derivadas dx(1), dx(2), ...
    % x0: vetor das condições iniciais do sistema
    
    fid = fopen('process.m', 'w');
    
    fprintf(fid, 'function sys = process(t,x,flag,u,p)\n');
    fprintf(fid, 'switch flag\n');
    
    % Caso 0: dinâmica
    fprintf(fid, '    case 0\n');
    fprintf(fid, '        %% Equações dinâmicas\n');
    
    % Extrair o número de estados
    n = length(dx_strings);
    
    % Definir as variáveis x1, x2,... e u1, u2,...
    for i=1:n
        fprintf(fid, '        x%d = x(%d);\n', i, i);
    end
    fprintf(fid, '        u1 = u(1);\n\n');
    
    fprintf(fid, '        %% Cálculo das derivadas\n');
    fprintf(fid, '        dx = zeros(%d,1);\n', n);
    for i=1:n
        % Trocar x(...) e u(...) para x1, x2,... e u1, u2,... na expressão
        expr = dx_strings{i};
        expr = regexprep(expr, 'x\((\d+)\)', 'x$1');
        expr = regexprep(expr, 'u\((\d+)\)', 'u$1');
        
        fprintf(fid, '        dx(%d) = %s;\n', i, expr);
    end
    
    fprintf(fid, '        sys = dx;\n\n');
    
    % Caso 5: condições iniciais
    fprintf(fid, '    case 5\n');
    fprintf(fid, '        sys = [%s];\n', num2str(x0));
    
    % Outros cases
    fprintf(fid, '    otherwise\n');
    fprintf(fid, '        sys = [];\n');
    fprintf(fid, 'end\n');
    
    fclose(fid);
    disp('Arquivo process.m criado com sucesso!');
end
