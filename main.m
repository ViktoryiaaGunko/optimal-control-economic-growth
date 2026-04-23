% ПАРАМЕТРЫ МОДЕЛИ 
T = 10;          
N = 30;          % Уменьшили количество точек для скорости fminsearch
dt = T / (N-1);  
t = linspace(0, T, N); 

alpha = 0.3;     
mu = 0.03;       
rho = 0.05;      
epsilon = 0.01;  
z0 = 2.0;        

% НАСТРОЙКИ ОПТИМИЗАЦИИ (MATLAB - fminsearch) 
% Чтобы ограничить u от 0 до 0.99, мы оптимизируем переменную v
% Связь: u = 0.99 / (1 + exp(-v))
u_init = 0.11; % Стартуем с магистрали
v0 = -log((1 - epsilon)/u_init - 1) * ones(1, N); % Обратное преобразование

% Настройки решателя
options = optimset('Display', 'iter', 'MaxIter', 5000, 'MaxFunEvals', 20000);

% Запуск оптимизации
v_opt = fminsearch(@(v) cost_function_unconstrained(v, N, dt, t, alpha, rho, mu, z0, epsilon), v0, options);

% ВОССТАНОВЛЕНИЕ РЕЗУЛЬТАТОВ 
% Возвращаем u_opt в диапазон [0, 0.99]
u_opt = (1 - epsilon) ./ (1 + exp(-v_opt));

% Прямой расчет z_opt по найденному u_opt
z_opt = zeros(1, N);
z_opt(1) = z0;
for i = 1:N-1
    z_opt(i+1) = z_opt(i) + dt * (u_opt(i) * z_opt(i)^alpha - mu * z_opt(i));
end

% ПОСТРОЕНИЕ ГРАФИКОВ 
figure('Name', 'Результаты оптимального управления');

% 1. График траектории капитала (z)
subplot(1,3,1);
plot(t, z_opt, 'b-', 'LineWidth', 2);
hold on;
yline(6.46, 'r--', 'Магистраль (6.46)', 'LineWidth', 1.5);
grid on;
title('Траектория капитала z(t)');
xlabel('Время, t');
ylabel('Капитал, z');

% 2. График управления (u)
subplot(1,3,2);
plot(t, u_opt, 'g-', 'LineWidth', 2);
hold on;
yline(0.11, 'r--', 'Магистраль (0.11)', 'LineWidth', 1.5);
ylim();
grid on;
title('Управление u(t)');
xlabel('Время, t');
ylabel('Доля инвестиций, u');

% 3. Фазовый портрет (z vs u)
subplot(1,3,3);
plot(z_opt, u_opt, 'k-', 'LineWidth', 2);
hold on;
plot(z_opt(1), u_opt(1), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'b');
plot(z_opt(end), u_opt(end), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
plot(6.46, 0.11, 'm*', 'MarkerSize', 10);
grid on;
title('Фазовый портрет');
xlabel('Капитал, z');
ylabel('Управление, u');
legend('Траектория', 'Начало (t=0)', 'Конец (t=T)', 'Стационарная точка');

