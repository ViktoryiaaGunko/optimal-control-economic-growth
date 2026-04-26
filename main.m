clear; clc; close all;

% параметры
alpha = 0.3;
mu = 0.03;
rho = 0.05;
epsilon = 0.01;
z0 = 2;              % Начальный капитал
N = 100;             % Количество точек

% Расчет магистрального решения
z_bar = (alpha / (rho + mu))^(1 / (1 - alpha));
u_bar = mu * z_bar^(1 - alpha);

T_list = [20, 30, 50];
colors = {'b', 'g', 'r'}; % Цвета для T=20, 30, 50

% графики
fig_z = figure('Name', 'Траектория капитала z(t)', 'Position', [100, 100, 800, 500]);
hold on; grid on;
title('Траектория капитала z(t)');
xlabel('Время t'); ylabel('Капитал z(t)');

fig_u = figure('Name', 'Управление u(t)', 'Position', [150, 150, 800, 500]);
hold on; grid on;
title('Оптимальное управление u(t)');
xlabel('Время t'); ylabel('Доля инвестиций u(t)');
ylim([0 0.3]);

fig_phase = figure('Name', 'Фазовый портрет', 'Position', [200, 200, 800, 500]);
hold on; grid on;
title('Фазовый портрет');
xlabel('Капитал z(t)'); ylabel('Скорость изменения dz/dt');

% цикл по горизонтам планирования
for k = 1:length(T_list)
    T = T_list(k);
    dt = T / (N - 1);
    t = linspace(0, T, N);

    % Начальное приближение и границы
    u_init = u_bar * ones(1, N);
    lb = zeros(1, N);
    ub = (1 - epsilon) * ones(1, N);

    % Настройки fmincon
    options = optimoptions('fmincon', ...
        'Algorithm', 'sqp', ...
        'Display', 'iter', ...
        'MaxFunctionEvaluations', 100000, ...
        'MaxIterations', 1000);

    % Запуск оптимизации
    fprintf('\n---> Решение задачи для T = %d...\n', T);
    u_opt = fmincon(@(u) cost_function(u, t, dt, z0, alpha, mu, rho), ...
                    u_init, [], [], [], [], lb, ub, [], options);

    % Восстановление траектории капитала z(t)
    z_opt = zeros(1, N);
    z_opt(1) = z0;
    for i = 1:N-1
        z_opt(i+1) = z_opt(i) + dt * (u_opt(i) * z_opt(i)^alpha - mu * z_opt(i));
    end

    % Вычисление производной dz/dt для фазового портрета
    zdot_opt = u_opt .* z_opt.^alpha - mu .* z_opt;
    
    % График z(t)
    figure(fig_z);
    plot(t, z_opt, 'Color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('T = %d', T));

    % График u(t)
    figure(fig_u);
    plot(t, u_opt, 'Color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('T = %d', T));

    % Фазовый портрет
    figure(fig_phase);
    plot(z_opt, zdot_opt, 'Color', colors{k}, 'LineWidth', 2, 'DisplayName', sprintf('T = %d', T));
end

% добавление магистралей и легенд 
% окно z(t)
figure(fig_z);
yline(z_bar, '--k', 'Магистраль \barz', 'LineWidth', 1.5, 'HandleVisibility', 'off', 'LabelHorizontalAlignment', 'left');
legend('Location', 'best');

% окно u(t)
figure(fig_u);
yline(u_bar, '--k', 'Магистраль \baru', 'LineWidth', 1.5, 'HandleVisibility', 'off', 'LabelHorizontalAlignment', 'left');
legend('Location', 'best');

% фазовое окно
figure(fig_phase);
plot(z_bar, 0, 'k*', 'MarkerSize', 10, 'LineWidth', 1.5, 'DisplayName', 'Стационарная точка');
legend('Location', 'best');

fprintf('\nВсе вычисления завершены. Графики построены в отдельных окнах!\n');

% функция
function J = cost_function(u, t, dt, z0, alpha, mu, rho)
    N = length(u);
    z = zeros(1, N);
    z(1) = z0;
    J = 0;

    for i = 1:N-1
        z(i+1) = z(i) + dt * (u(i) * z(i)^alpha - mu * z(i));
        
        if z(i) <= 0 || u(i) >= 1
            J = 1e10; 
            return;
        end
        
        utility = exp(-rho * t(i)) * (log(1 - u(i)) + alpha * log(z(i)));
        J = J + utility * dt; 
    end
    
    if z(end) <= 0 || u(end) >= 1
        J = 1e10;
        return;
    end
    
    utility_end = exp(-rho * t(end)) * (log(1 - u(end)) + alpha * log(z(end)));
    J = J + utility_end * dt;
    
    J = -J; 
end
