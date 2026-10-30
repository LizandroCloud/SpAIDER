function u = get_u(t, ti, ui)
    % Assume: length(ti) == length(ui)
    n = length(ti);
    
    for i = 1:n
        if i < n
            if t < ti(i)
                u = ui(i);
                return;
            end
        else
            % Último intervalo: t >= ti(end)
            u = ui(end);
            return;
        end
    end
end
