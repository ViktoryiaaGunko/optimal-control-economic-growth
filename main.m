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

% --- НАСТРОЙКИ ОПТИМИЗАЦИИ (MATLAB - fminsearch) ---
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
