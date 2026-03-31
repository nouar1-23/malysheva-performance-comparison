function [d1, d2, d3,d4] = con_d(a, b, n)
%  Сравнительный анализ времени работы трёх вариантов метода Малышева.
% 
% Входные параметры:
%   a - Начальное значение параметра λ (радиус)
%   b - Конечное значение параметра λ
%   n - Размерность случайной матрицы A
%
% Выходные параметры:
%   d1, d2, d3 - Время выполнения (в секундах) для каждого алгоритма

    % --- Инициализация общих параметров ---
    r1_val = a;          % Текущее значение радиуса
    h      = 0.1;        % Шаг по параметру
    w_max  = 1e10;       % Ограничение нормы решения
    u_max  = 1e10;       % Ограничение обусловленности
    ip     = 1e-13;      % Точность вычислений
    
    A      = rand(n, n); % Генерируем случайную матрицу
    I      = eye(n);     % Единичная матрица
    
    % Расчет максимально необходимого числа итераций m0
    m0 = round(log2(- (1 + w_max) * log(ip / ((2 + 2 * ip) * sqrt(w_max)))));

    % ================================================================
   % % АЛГОРИТМ 1: Фиксированное число итераций без проверки проектора
    % ================================================================
    i = 1; 
    tic; % Старт таймера 1
    while (r1_val <= b)
        skip_step = false; 
        A0 = A' / r1_val;
        B0 = I;
        d  = 0;
        
        while d <= m0 / 3
            d = d + 1;
            for i1 = 1:3
                S_temp = [-B0; A0];
                [Q_mat, ~] = qr(S_temp);
                QQ     = Q_mat';
                A0     = QQ(n+1:2*n, 1:n) * A0;
                B0     = QQ(n+1:2*n, n+1:2*n) * B0;
            end
            if (cond(A0 - B0) > u_max)
                skip_step = true; break; 
            end
        end
        
        if ~skip_step
            inv_su = inv(B0 + A0);
            H      = inv_su * inv_su';
            m1(i)  = max(abs(H), [], 'all');
            e1(i)  = r1_val;
            i      = i + 1;
        end
        r1_val = r1_val + h;
    end
    d1 = toc; % Финиш таймера 1

    % ================================================================
    %% АЛГОРИТМ 2: С условием сходимости по свойству проектора (p^2 = p)
    % ================================================================
    i = 1; 
    r1_val = a;
    tic; % Старт таймера 2
    while (r1_val <= b)
        skip_step = false; 
        A0 = A' / r1_val;
        B0 = I;
        p  = 2 * I; % Инициализация проектора
        d  = 0;
        
        while (max(abs(p^2 - p), [], 'all') > ip  && d <= m0/3)
            d = d + 1;
            for i1 = 1:3
                S_temp = [-B0; A0];
                [Q_mat, ~] = qr(S_temp);
                QQ     = Q_mat';
                A0     = QQ(n+1:2*n, 1:n) * A0;
                B0     = QQ(n+1:2*n, n+1:2*n) * B0;
            end
            if cond(A0 - B0) > u_max
                skip_step = true; break; 
            end
            p = -inv(A0 - B0) * B0; 
        end
        
        if ~skip_step && ~isempty(p) && ~any(isnan(p(:)))
            inv_su = inv(B0 + A0);
            H      = inv_su * inv_su';
            m2(i)  = max(abs(H), [], 'all');
            e2(i)  = r1_val;
            i      = i + 1;
        end
        r1_val = r1_val + h;
    end
    d2 = toc;

    % ================================================================
    %% АЛГОРИТМ 3: С условием сходимости по изменению нормы решения H
    % ================================================================
    i = 1; 
    r1_val = a;
    tic; % Старт таймера 3
    while (r1_val <= b)
        skip_step = false; 
        A0 = A' / r1_val;
        B0 = I;
        H = zeros(n);
        HH = I; 
        d  = 0;
        
        while (max(abs(H - HH), [], 'all') > ip * max(abs(H), [], 'all') && d <= m0/3)
            d = d + 1;
            for i1 = 1:3
                S_temp = [-B0; A0];
                [Q_mat, ~] = qr(S_temp);
                QQ     = Q_mat';
                A0     = QQ(n+1:2*n, 1:n) * A0;
                B0     = QQ(n+1:2*n, n+1:2*n) * B0;
            end
            if (cond(A0 - B0) > u_max), skip_step = true; break; end
            
            % Текущее решение
            inv_su = inv(B0 + A0);
            H = inv_su * inv_su';
            
            % Прогноз следующего шага для проверки сходимости
            S_temp = [-B0; A0];
            [Q_mat, ~] = qr(S_temp);
            QQ     = Q_mat';
            A1     = QQ(n+1:2*n, 1:n) * A0;
            B1     = QQ(n+1:2*n, n+1:2*n) * B0;
            inv_su_next = inv(A1 + B1);
            HH = inv_su_next * inv_su_next';
        end
        
        if ~skip_step
            m3(i) = max(abs(HH), [], 'all');
            e3(i) = r1_val;
            i     = i + 1;
        end
        r1_val = r1_val + h;
    end
    d3 = toc;
    % ================================================================
    %% АЛГОРИТМ 4: С условием сходимости  проектора ||p_j - p_{j-1}||
    % ================================================================
    i = 1; 
    r1_val = a;
    tic; % Старт таймера 2
    while (r1_val <= b)
        skip_step = false; 
        A0 = A' / r1_val;
        B0 = I;
        p1  = 2 * I; % Инициализация проектора
        p2=I;
        d  = 0;
        
        while (max(abs(p2 - p1), [], 'all') > ip  && d <= m0/3)
            d = d + 1;
            for i1 = 1:2
                S_temp = [-B0; A0];
                [Q_mat, ~] = qr(S_temp);
                QQ     = Q_mat';
                A0     = QQ(n+1:2*n, 1:n) * A0;
                B0     = QQ(n+1:2*n, n+1:2*n) * B0;
            end
            if cond(A0 - B0) > u_max
                skip_step = true; break; 
            end
            p1 = -inv(A0 - B0) * B0; 
            S_temp = [-B0; A0];
                [Q_mat, ~] = qr(S_temp);
                QQ     = Q_mat';
                A0     = QQ(n+1:2*n, 1:n) * A0;
                B0     = QQ(n+1:2*n, n+1:2*n) * B0;
                p2 = -inv(A0 - B0) * B0;
        end
        
        if ~skip_step && ~isempty(p2) && ~any(isnan(p2(:)))
            inv_su = inv(B0 + A0);
            H      = inv_su * inv_su';
            m4(i)  = max(abs(H), [], 'all');
            e4(i)  = r1_val;
            i      = i + 1;
        end
        r1_val = r1_val + h;
    end
    d4 = toc;

    %% --- Визуализация результатов ---
    figure;
    plot(e1, log10(m1), 'blue', 'LineWidth', 5); hold on;
    plot(e2, log10(m2), 'red', 'LineWidth', 3);
    plot(e3, log10(m3), 'k', 'LineWidth', 1);
    
    xlabel('|\lambda|');
    ylabel('log_{10}(||H||)');
    title(['Сравнение времени: Алг1=', num2str(d1), 's, Алг2=', num2str(d2), 's, Алг3=', num2str(d3), 's']);
    legend('Алгоритм 1', 'Алгоритм 2', 'Алгоритм 3');
    grid on;
end
