clear;
m_2 = 0.1; 
m_1 = 0.1; 

k_1 = 17/15;
b_2 = 17/15;

k_2 = 73/100;
b_1 = 73/100;

k_3 = 2/3 * k_1;
s = tf("s");
sysc = k_3 / ((m_2*s^2 + k_2 + k_3 + b_2*s)*(m_1*s^2+k_1+b_1*s+k_3) - k_3^2);

fb = bandwidth(sysc);
Ts= 0.05*2*pi/fb;
sysd = c2d(sysc,Ts,'zoh');

%% ELS

main_title = "ELS";
num_samples = 10000;
noise_variance = 0.01;

% step_mag = 1;
% t = 0:num_samples-1;
% u = step_mag * t>=0;

input_variance = 5;
u = sqrt(input_variance) * randn(1, num_samples);

[num, den]=tfdata(sysd,'v');

P_zero = 20000* eye(10);
theta_hat_zero = ones(10,1)*0;
noise_poly = [0.5,0.25];
theta_real = [den(2:end), num(2:end)].';



rng(50)
clear y_real y_predicted theta_hat phi phi_t P_new K_t
y_real = zeros(num_samples,1);
y_predicted = zeros(num_samples,1);
theta_hat = zeros(4+4+2, num_samples); %% depends on dynamics
phi = zeros(num_samples, 4+4+3);

epsilon = zeros(num_samples,1);

[num, den]=tfdata(sysd,'v');

theta_real = [den(2:end), num(2:end), noise_poly].';
noise_variance = noise_variance; % system noise variance
noise = sqrt(noise_variance) * randn(1, num_samples);

P_old = P_zero;
theta_hat_zero = theta_hat_zero;

for i=1:num_samples
    % generate output and phi_t
    switch i
    case 1 % at time t=1
        phi_t_without_noise = [zeros(1,4), zeros(1,4)].';
        noise_t = zeros(1,3);
        epsilon_t = zeros(1,2);
    case 2
        phi_t_without_noise = [-y_real(i-1),zeros(1,3), u(i-1),zeros(1,3)].';
        noise_t = [noise(i), noise(i-1),zeros(1,1)];
        epsilon_t = [epsilon(i-1),zeros(1,1)];
    case 3
        phi_t_without_noise = [-y_real(i-1:-1:i-2).',zeros(1,2), u(i-1:-1:i-2),zeros(1,2)].';
        noise_t = [noise(i), noise(i-1),noise(i-2)];
        epsilon_t = [epsilon(i-1),epsilon(i-2)];
    case 4
        phi_t_without_noise = [-y_real(i-1:-1:i-3).',zeros(1,1), u(i-1:-1:i-3),zeros(1,1)].';
        noise_t = [noise(i), noise(i-1),noise(i-2)];
        epsilon_t = [epsilon(i-1),epsilon(i-2)];
    otherwise
        phi_t_without_noise = [-y_real(i-1:-1:i-4).' u(i-1:-1:i-4)].';
        noise_t = [noise(i), noise(i-1),noise(i-2)];
        epsilon_t = [epsilon(i-1),epsilon(i-2)];
    end
    
    y_real(i) = (phi_t_without_noise.' * [den(2:end), num(2:end)].') + noise_t * [1; theta_real(9:10)];
    
    % approximating the new noise

    phi_t = [phi_t_without_noise; epsilon_t.'];
    P_new = P_old - (P_old * (phi_t * phi_t.') * P_old)/(1+phi_t.'*P_old*phi_t);
    K_t = (P_old * phi_t)/(1+phi_t.'*P_old*phi_t);
    

    switch i
        case 1 % at time t=1
            theta_hat(:, i) = theta_hat_zero + K_t * (y_real(i) - phi_t.' * theta_hat_zero) ;
            epsilon(i) = 0; 
        otherwise
            theta_hat(:, i) = theta_hat(:, i-1) + K_t * (y_real(i) - phi_t.' * theta_hat(:, i-1));
            % updating epsilon
            epsilon(i) = y_real(i) - phi_t.' * theta_hat(:, i);
            
    end
    y_predicted(i) = phi_t.' * theta_hat(:, i);
    
    P_old = P_new;
end

error_y = y_predicted-y_real;
SSE_y = norm(error_y,2)^2;
mse_y = immse(y_real, y_predicted);

error_theta = theta_real * ones(1,num_samples) - theta_hat;
SSE_theta = zeros(8,1);
mse_theta = zeros(8,1);
for k = 1:length(theta_real)
    SSE_theta(k) = norm(error_theta(k,:),2)^2;
    mse_theta(k) = mean(error_theta(k,:).^2);
end


% figures
f1 = figure();
f1.Position = [-1000 0 1000 500];

plot(1:num_samples, y_real, "--b", 'DisplayName','real output');
xlabel("sample time");
hold on;
plot(1:num_samples, y_predicted, 'DisplayName','predicted output');
legend('Location','best');
title("Output")
saveas(gcf,'images/q2/' + main_title+ "_summary" + '.jpeg')
% close all



f2 = figure();
f2.Position = [-1500 -200 1000 1300];

for i = 1:length(theta_hat(:,1))
    title_text = "θ_%d";
    subplot(5,2,i);
    plot(1:num_samples, theta_hat(i,:), 'DisplayName','Predicted')
    hold on;
    plot(1:num_samples, ones(size(1:num_samples)) * theta_real(i) , 'DisplayName','Real')
    title(sprintf(title_text, i));
    legend('Location','best');
    xlabel("sample number")
end

saveas(gcf,'images/q2/' + main_title + '.jpeg')
% close all

f3 = figure();
f3.Position = [-1000 0 1000 500];
sysd_predicted = tf([0, theta_hat(5:8, num_samples).'], [1 , theta_hat(1:4, num_samples).'], sysd.Ts); 
bode(sysd_predicted);
hold on; 
h = bodeplot(sysd);
setoptions(h,'FreqUnits','Hz','PhaseVisible','off');
legend("Predicted system", "Real system");
saveas(gcf,'images/q2/' + main_title+ "_bode" + '.jpeg')
% close all

% [SSE_y, SSE_theta, mse_y, mse_theta, theta_hat] = ELS(u, num_samples, sysd,noise_variance, theta_hat_zero, P_zero, noise_poly, main_title);
writematrix(round(theta_hat(:,num_samples),6),"theta_hat.txt");
writematrix(round([SSE_y;mse_y],6), "SSE_y.txt");
writematrix(round([SSE_theta;mse_theta],6), "SSE_theta.txt");
% 
