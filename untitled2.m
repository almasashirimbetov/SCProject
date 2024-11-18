%% Parameters and Initial Setup
% Define constants
m = 1; % Mass (kg)
k = 4; % Spring constant (N/m)

% Natural frequency
omega_n = sqrt(k / m);

% Time span
t_span = [0, 10]; % Simulation time span

%% Analytical and Numerical Solution for Default Case
% Default parameters
zeta_default = 0.2; % Default damping ratio (underdamped)
x0_default = 1;     % Initial displacement (m)
v0_default = 0;     % Initial velocity (m/s)

% Analytical solution for default case (underdamped)
omega_d = omega_n * sqrt(1 - zeta_default^2); % Damped natural frequency
x_analytical_default = @(t) exp(-zeta_default * omega_n * t) .* ...
    (x0_default * cos(omega_d * t) + ...
    (v0_default + zeta_default * omega_n * x0_default) / omega_d * sin(omega_d * t));

% Numerical solution for default case
damped_oscillator = @(t, y) [y(2); -2 * zeta_default * omega_n * y(2) - omega_n^2 * y(1)];
[t_num_default, y_num_default] = ode45(damped_oscillator, t_span, [x0_default; v0_default]);

% Plot the default case
figure;
plot(t_num_default, y_num_default(:, 1), 'b-', 'LineWidth', 1.5); hold on;
fplot(x_analytical_default, t_span, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Displacement (m)');
legend('Numerical Solution', 'Analytical Solution');
title('Default Case: Underdamped (\zeta = 0.2)');
grid on;

% Display max error for default case
x_analytical_values_default = arrayfun(x_analytical_default, t_num_default);
default_error = abs(x_analytical_values_default - y_num_default(:, 1));
fprintf('Default Case Max Error: %.5f\n', max(default_error));

%% Varying Initial Conditions
% Damping ratios for testing
zeta_values = [0.2, 1.0, 2.0]; % Underdamped, critically damped, overdamped

% Initial conditions for testing
x0_values = [1, 0.5, -1]; % Initial displacements
v0_values = [0, 1, -0.5]; % Initial velocities

% Loop over damping ratios
for zeta = zeta_values
    fprintf('\nTesting for zeta = %.1f\n', zeta);
    
    % Loop over initial conditions
    for x0 = x0_values
        for v0 = v0_values
            fprintf('  Initial conditions: x0 = %.1f, v0 = %.1f\n', x0, v0);
            
            % Define the ODE for the current test case
            damped_oscillator = @(t, y) [y(2); -2 * zeta * omega_n * y(2) - omega_n^2 * y(1)];
            y0 = [x0; v0];
            
            % Solve numerically with ode45
            [t_num, y_num] = ode45(damped_oscillator, t_span, y0);
            x_numerical = y_num(:, 1); % Numerical displacement
            
            % Define the analytical solution for the current test case
            if zeta < 1
                % Underdamped case
                omega_d = omega_n * sqrt(1 - zeta^2);
                x_analytical = @(t) exp(-zeta * omega_n * t) .* ...
                    (x0 * cos(omega_d * t) + ...
                    (v0 + zeta * omega_n * x0) / omega_d * sin(omega_d * t));
            elseif zeta == 1
                % Critically damped case
                x_analytical = @(t) exp(-omega_n * t) .* (x0 + (v0 + omega_n * x0) .* t);
            else
                % Overdamped case
                lambda1 = -zeta * omega_n + omega_n * sqrt(zeta^2 - 1);
                lambda2 = -zeta * omega_n - omega_n * sqrt(zeta^2 - 1);
                C1 = (v0 - lambda2 * x0) / (lambda1 - lambda2);
                C2 = x0 - C1;
                x_analytical = @(t) C1 * exp(lambda1 * t) + C2 * exp(lambda2 * t);
            end
            
            % Evaluate the analytical solution at the numerical time points
            x_analytical_values = arrayfun(x_analytical, t_num);
            
            % Plot numerical vs. analytical results
            figure;
            plot(t_num, x_numerical, 'b-', 'LineWidth', 1.5); hold on;
            plot(t_num, x_analytical_values, 'r--', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Displacement (m)');
            legend('Numerical Solution', 'Analytical Solution');
            title(sprintf('Damped Oscillator (zeta = %.1f, x0 = %.1f, v0 = %.1f)', zeta, x0, v0));
            grid on;
            
            % Calculate and display max error for the current test case
            error = abs(x_analytical_values - x_numerical);
            fprintf('    Max Error: %.5f\n', max(error));
            
            % Plot the error over time
            figure;
            plot(t_num, error, 'k-', 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Error (m)');
            title(sprintf('Error (zeta = %.1f, x0 = %.1f, v0 = %.1f)', zeta, x0, v0));
            grid on;
        end
    end
end