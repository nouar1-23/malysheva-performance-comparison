function [] = tim(n)

%Размеры матриц, для которых требуется вычислить время, рассматриваются на интервале от 10 до n с шагом 10.
    k = 10:10:n;
    m = length(k);
c=15;
    for j = 1:m
        ww1 = []; ww2 = []; ww3 = []; ww4 = [];
        % Вычисление времени, необходимого для трёх случайных матриц одного и того же размера.
        for x = 1:c
            [w1, w2, w3,w4] = con_d_tim(1, 2, k(j));
            ww1(x) = w1;
            ww2(x) = w2;
            ww3(x) = w3;
            ww4(x) = w4;
        end
        % Вычисление среднего арифметического времени, затраченного в трёх повторных запусках.
        z1{j} = ww1;
        z2{j} = ww2;
        z3{j} = ww3;
        z4{j} = ww4;
    end
% Построение графиков, отображающих время, необходимое для каждого варианта.
    figure; hold on;

for j = 1:length(k)
        plot(ones(1,c)*k(j), z1{j}, 'b+', 'MarkerSize', 6, 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), z4{j}, 'k*', 'MarkerSize', 5, 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), z2{j}, 'ro', 'MarkerSize', 4, 'HandleVisibility', 'off');
        plot(ones(1,c)*k(j), z3{j}, 'g*', 'MarkerSize', 4, 'HandleVisibility', 'off');

end
    mean_z1 = cellfun(@mean, z1);
    mean_z2 = cellfun(@mean, z2);
    mean_z3 = cellfun(@mean, z3);
     mean_z4 = cellfun(@mean, z4);

    plot(k, mean_z1, 'b', 'LineWidth', 1.4); hold on;
    plot(k, mean_z4, 'k', 'LineWidth', 1.4);
    plot(k, mean_z2, 'r', 'LineWidth', 1.4);
    plot(k, mean_z3, 'g', 'LineWidth', 1.4);
    

    xlabel('Размерность матрицы, n');
    ylabel('Время (секунды)');
    set(gca, 'FontSize', 25);
    legend('Алг.1','Алг.2, Вар.1','Алг.2, Вар.2','Алг.2, Вар.3','FontSize', 25, 'Location', 'northwest');
% figure
% hold on
%     plot(k, z2, 'r', 'LineWidth', 1.4);
%     plot(k, z3, 'g', 'LineWidth', 1.4);
%     plot(k, z4, 'k', 'LineWidth', 1.4);
% 
%     xlabel('Размерность матрицы, n');
%     ylabel('Время (секунды)');
%     set(gca, 'FontSize', 25);
%      legend('Алг.2, Вар.1','Алг.2, Вар.2','Алг.2, Вар.3','FontSize', 25, 'Location', 'northwest');
