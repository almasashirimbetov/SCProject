% Parameters
m = 100;       % Mass (kg)
k = 400;       % Spring constant (N/m)
c = 0.2;     % Damping coefficient (N·s/m)
x0 = 1;      % Initial displacement (m)
v0 = 0;      % Initial velocity (m/s)

% Derived quantities
omega_n = sqrt(k/m);       % Natural frequency
zeta = c / (2 * sqrt(m*k)); % Damping ratio

% Symbolic solution
syms t
omega_d = omega_n * sqrt(1 - zeta^2); % Damped natural frequency

if zeta < 1
    % Underdamped case
    x_analytical = exp(-zeta * omega_n * t) * (x0 * cos(omega_d * t) + ...
                  (v0 + zeta * omega_n * x0) / omega_d * sin(omega_d * t));
elseif zeta == 1
    % Critically damped case
    x_analytical = exp(-omega_n * t) * (x0 + (v0 + omega_n * x0) * t);
else
    % Overdamped case
    lambda1 = -zeta * omega_n + omega_n * sqrt(zeta^2 - 1);
    lambda2 = -zeta * omega_n - omega_n * sqrt(zeta^2 - 1);
    C1 = (v0 - lambda2 * x0) / (lambda1 - lambda2);
    C2 = x0 - C1;
    x_analytical = C1 * exp(lambda1 * t) + C2 * exp(lambda2 * t);
end

% Convert symbolic solution to function handle
x_analytical_fn = matlabFunction(x_analytical, 'Vars', t);

% Define the ODE
damped_oscillator = @(t, y) [y(2); -2 * zeta * omega_n * y(2) - omega_n^2 * y(1)];

% Time span and initial conditions
t_span = [0, 10];
y0 = [x0; v0];

% Solve ODE numerically
[t_num, y_num] = ode45(damped_oscillator, t_span, y0);

% Extract displacement
x_numerical = y_num(:, 1);

% Evaluate analytical solution at numerical time points
x_analytical_values = arrayfun(x_analytical_fn, t_num);

% Plot both solutions
figure;
plot(t_num, x_numerical, 'b-', 'LineWidth', 1.5); hold on;
plot(t_num, x_analytical_values, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Displacement (m)');
legend('Numerical Solution', 'Analytical Solution');
title('Damped Harmonic Oscillator: Analytical vs Numerical');
grid on;

% Compute error
error = abs(x_analytical_values - x_numerical);

% Plot error
figure;
plot(t_num, error, 'k-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Error (m)');
title('Error Between Analytical and Numerical Solutions');
grid on;

% Example 1: Underdamped
c = 0.2; % Damping coefficient
run_simulation(m, c, k, x0, v0);

% Example 2: Critically damped
c = 2 * sqrt(m * k);
run_simulation(m, c, k, x0, v0);

% Example 3: Overdamped
c = 8;
run_simulation(m, c, k, x0, v0);

function run_simulation(m, c, k, x0, v0)
    omega_n = sqrt(k / m);
    zeta = c / (2 * sqrt(m * k));
    fprintf('Running simulation for zeta = %.2f\n', zeta);

    % Define ODE
    damped_oscillator = @(t, y) [y(2); -2 * zeta * omega_n * y(2) - omega_n^2 * y(1)];

    % Time span and initial conditions
    t_span = [0, 10];
    y0 = [x0; v0];

    % Solve numerically
    [t_num, y_num] = ode45(damped_oscillator, t_span, y0);
    x_numerical = y_num(:, 1);

    % Plot numerical solution
    figure;
    plot(t_num, x_numerical, 'b-', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Displacement (m)');
    title(sprintf('Damped Harmonic Oscillator (zeta = %.2f)', zeta));
    grid on;
end