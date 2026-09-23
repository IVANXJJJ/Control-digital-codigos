
% TAREA 1 DE CONTROL DIGITAL: MOTOR DC Y  CUANTIZACIÓN
% Autor: Iván Jara
% Planta: G(s) = Km / (s * (Tm*s + 1))  |  Km = 0.5, Tm = 0.1 s

clc; clear; close all;


%% PARTE 1: RESOLUCIÓN Y ERROR DE CUANTIZACIÓN

V_max = 5;
V_min = -5;
V_fsr = V_max - V_min; % Rango de Escala Completa = 10 V
bits_ADC = 10;
bits_DAC = 8;

% 1. Cálculos de resolución (LSB) y error máximo (+- LSB/2)
delta_ADC = V_fsr / (2^bits_ADC);        % 10 / 1024 = 9.7656 mV
delta_DAC = V_fsr / (2^bits_DAC);        % 10 / 256  = 39.0625 mV
err_max_ADC = delta_ADC / 2;             % +- 4.8828 mV
err_max_DAC = delta_DAC / 2;             % +- 19.5313 mV

fprintf('--------------------------------------------------------\n');
fprintf(' 1. RESOLUCIÓN DE CONVERTIDORES Y ERROR DE CUANTIZACIÓN\n');
fprintf('-----------------------------------------------------\n');
fprintf('- ADC (10 bits): Res = %.4f V (%.2f mV) | Err Máx = +-%0.4f V\n', delta_ADC, delta_ADC*1000, err_max_ADC);
fprintf('- DAC ( 8 bits): Res = %.4f V (%.2f mV) | Err Máx = +-%0.4f V\n\n', delta_DAC, delta_DAC*1000, err_max_DAC);

fprintf('INTERPRETACIÓN FÍSICA PARA EL LAZO DE CONTROL:\n');
fprintf('1. Lectura (ADC): Cambios en la velocidad menores a %.2f mV no generan cambio\n', delta_ADC*1000);
fprintf('   en el código digital. Representa la zona muerta de medición del sensor.\n');
fprintf('2. Acción (DAC): El microcontrolador solo puede modificar la tensión del motor en\n');
fprintf('   pasos discretos de %.2f mV, lo que impide correcciones infinitesimales.\n\n', delta_DAC*1000);


%% PARTE 2: DISEÑO DEL PID DIGITAL 

Km = 0.5;      % Ganancia del motor [rad/(V*s)]
Tm = 0.1;      % Constante de tiempo del motor [s]
T  = 0.01;     % Período de muestreo [s] (10 ms)

fprintf('----------------------------------------------------\n');
fprintf(' 2. DISEÑO Y JUSTIFICACIÓN DEL CONTROLADOR PID\n');
fprintf('--------------------------------------------------------\n');
fprintf('JUSTIFICACIÓN DE T = 10 ms:\n');
fprintf('La constante del motor es Tm = 100 ms. La relación Tm/T = 100ms/10ms = 10,\n');
fprintf('lo cual cumple el criterio de diseño discreto (10 <= Tm/T <= 20).\n\n');

% 1. Requisitos de especificación
Mp_des = 0.05; % Sobrepaso deseado (<= 5%)
ts_des = 0.50; % Tiempo de asentamiento deseado (< 500 ms)

% 2. Cálculo analítico base de 2do orden
zeta_base = sqrt(log(Mp_des)^2 / (pi^2 + log(Mp_des)^2)); % 0.6901
wn_base   = 4 / (zeta_base * ts_des);                     % 11.59 rad/s

% 3. Criterio de Margen por Ceros del PID y Retardo ZOH
gamma_robustez = 1.20; 
zeta = zeta_base * gamma_robustez; % zeta = 0.8281 (~0.83)
wn   = ceil(wn_base);               % 12 rad/s

% 4. Matriz de ganancias calculadas por igualación de coeficientes
Kp = (wn^2 * Tm) / Km;                   % Kp = 28.80
Kd = ((2 * zeta * wn * Tm) - 1) / Km;     % Kd = 1.987
Ki = 0.50;                                % Término integral bajo para asegurar e_ss = 0

fprintf('Parámetros base teórico: zeta_base = %.4f | wn_base = %.2f rad/s\n', zeta_base, wn_base);
fprintf('Parámetros con margen  : zeta = %.4f (gamma=%.2f) | wn = %.2f rad/s\n', zeta, gamma_robustez, wn);
fprintf('Ganancias diseñadas    -> Kp = %.2f | Ki = %.2f | Kd = %.2f\n\n', Kp, Ki, Kd);

% Modelo discreto exacto usando ZOH
s = tf('s');
Gp = Km / (s * (Tm * s + 1));
Gz = c2d(Gp, T, 'zoh');
Cz = pid(Kp, Ki, Kd, 'Ts', T, 'IFormula', 'Backward', 'DFormula', 'Backward');
Tz_ideal = feedback(Cz * Gz, 1);


%% PARTE 3: VERIFICACIÓN DE ESPECIFICACIONES IDEALES

t_sim = 0:T:1.5;
N = length(t_sim);
r = ones(size(t_sim)); 
[y_ideal, t_out] = step(Tz_ideal, t_sim);
info_ideal = stepinfo(y_ideal, t_out, 1);

fprintf('-----------------------------------------------------\n');
fprintf(' 3. VERIFICACIÓN DE MÉTRICAS EN LAZO CERRADO (IDEAL)\n');
fprintf('----------------------------------------------------------\n');
fprintf('  Sobrepaso Máximo (Mp)  : %.2f %%  (Requisito <= 5%%)  -> ', info_ideal.Overshoot);
if info_ideal.Overshoot <= 5, fprintf('CUMPLE\n'); else, fprintf('NO CUMPLE\n'); end
fprintf('  Tiempo Asentamiento(ts): %.2f ms (Requisito < 500 ms) -> ', info_ideal.SettlingTime*1000);
if (info_ideal.SettlingTime*1000) < 500, fprintf('CUMPLE\n\n'); else, fprintf('NO CUMPLE\n\n'); end


%% PARTE 4: SIMULACIÓN CON CUANTIZACIÓN Y MITIGACIÓN 

[num_z, den_z] = tfdata(Gz, 'v');

% --- Caso A: PID Estándar (Sin filtro) ---
y_q1 = zeros(1, N); u_q1 = zeros(1, N); e_q1 = zeros(1, N);
integral_e1 = 0;

% --- Caso B: PID Mitigado (Con Filtro Derivativo) ---
y_q2 = zeros(1, N); u_q2 = zeros(1, N); e_q2 = zeros(1, N);
integral_e2 = 0; D_filt = 0;
N_filter = 8; 
alpha = (Kd/N_filter) / ((Kd/N_filter) + T); 

for k = 3:N
    % *** CASO A: SIN FILTRO ***
    y_adc1 = round(y_q1(k-1) / delta_ADC) * delta_ADC;
    e_q1(k) = r(k) - y_adc1;
    P1 = Kp * e_q1(k);
    integral_e1 = integral_e1 + e_q1(k) * T;
    I1 = Ki * integral_e1;
    D1 = (Kd / T) * (e_q1(k) - e_q1(k-1)); 
    u_calc1 = P1 + I1 + D1;
    u_sat1 = max(min(u_calc1, V_max), V_min);
    u_q1(k) = round(u_sat1 / delta_DAC) * delta_DAC;
    y_q1(k) = -den_z(2)*y_q1(k-1) - den_z(3)*y_q1(k-2) + num_z(2)*u_q1(k-1) + num_z(3)*u_q1(k-2);
    
    % *** CASO B: CON FILTRO DERIVATIVO ***
    y_adc2 = round(y_q2(k-1) / delta_ADC) * delta_ADC;
    e_q2(k) = r(k) - y_adc2;
    P2 = Kp * e_q2(k);
    integral_e2 = integral_e2 + e_q2(k) * T;
    I2 = Ki * integral_e2;
    D_filt = alpha * D_filt + (1 - alpha) * (Kd / T) * (e_q2(k) - e_q2(k-1));
    u_calc2 = P2 + I2 + D_filt;
    u_sat2 = max(min(u_calc2, V_max), V_min);
    u_q2(k) = round(u_sat2 / delta_DAC) * delta_DAC;
    y_q2(k) = -den_z(2)*y_q2(k-1) - den_z(3)*y_q2(k-2) + num_z(2)*u_q2(k-1) + num_z(3)*u_q2(k-2);
end

fprintf('====================================================\n');
fprintf(' 4. ANÁLISIS DEL EFECTO DE CUANTIZACIÓN Y MITIGACIÓN\n');
fprintf('--------------------------------------------------------------\n');
fprintf('- Término más afectado: DERIVATIVO (Kd).\n');
fprintf('  Un salto mínimo de 1 LSB del ADC (%.2f mV) genera un pico derivativo\n', delta_ADC*1000);
fprintf('  de (Kd/T)*delta_ADC = (%.2f / 0.01) * %.4f = %.2f V en u[k].\n', Kd, delta_ADC, (Kd/T)*delta_ADC);
fprintf('  Esto satura el DAC continuamente y genera oscilaciones de alta frecuencia.\n');
fprintf('- Mitigación implementada: Filtro Pasa-Bajos en el término derivativo (N = %d).\n', N_filter);
fprintf('  Suaviza los escalones del ADC eliminando el chattering en la señal del DAC.\n\n');


%% GRAFICACIÓN DE RESULTADOS

figure('Name', 'Control Digital Motor DC', 'Color', [1 1 1]);

% Subplot 1: Respuesta del Sistema
subplot(2,1,1);
plot(t_sim, y_ideal, 'b-', 'LineWidth', 1.8, 'DisplayName', 'Ideal (Sin Cuantizar)'); hold on;
plot(t_sim, y_q1, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Cuantizado (Sin Filtro)');
plot(t_sim, y_q2, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Cuantizado (Con Filtro N=8)');
plot(t_sim, r, 'k:', 'LineWidth', 1.2, 'DisplayName', 'Referencia');

% Líneas de referencia para métricas de evaluación
yline(1.05, 'm--', 'Límite +5% Mp', 'LineWidth', 1.0, 'HandleVisibility', 'off');
xline(0.50, 'k--', 'Límite ts = 500 ms', 'LineWidth', 1.0, 'HandleVisibility', 'off');

title('Respuesta Temporal del Motor DC (Velocidad Angular)');
xlabel('Tiempo (s)'); ylabel('Velocidad (rad/s)');
legend('Location', 'SouthEast'); grid on;

% Subplot 2: Esfuerzo de Control u[k]
subplot(2,1,2);
stairs(t_sim, u_q1, 'r--', 'LineWidth', 1.0, 'DisplayName', 'DAC Sin Filtro (Chattering)'); hold on;
stairs(t_sim, u_q2, 'g-', 'LineWidth', 1.4, 'DisplayName', 'DAC Con Filtro (Estable)');
title('Señal de Control u[k] Generada por el DAC (8 bits)');
xlabel('Tiempo (s)'); ylabel('Voltaje (V)');
ylim([-1 6]);
legend('Location', 'NorthEast'); grid on;
