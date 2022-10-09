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

%figure();
%bode(sysd)
%hold on
%bode(sysc)


rng(50) % random seed
%% white noise as the input 


num_samples = 500;
input_variance = 1;
u = sqrt(input_variance) * randn(1, num_samples);
[mse, sysd_predicted, theta_hat] = least_sqr(u, num_samples, sysd);
writematrix(round(theta_hat,6));

%% step as the input signal

num_samples = 500;
step_mag = 1;
t = 0:num_samples-1;
u = step_mag * t>=0;
[mse, sysd_predicted, theta_hat] = least_sqr(u, num_samples, sysd);
writematrix(round(theta_hat,6));

%%  pulse as the input signal

num_samples = 500;
puls_mag = 1;
t = 0:num_samples-1;
u = zeros(1, num_samples);
u(1:9) = puls_mag;
[mse, sysd_predicted, theta_hat] = least_sqr(u, num_samples, sysd);
writematrix(round(theta_hat,6));

%% sin

num_samples = 500;
sinus_mag = 1;
siuns_frequency = 0.5; 
t = 0:num_samples-1;
u = sinus_mag * sin(siuns_frequency * t);
[mse, sysd_predicted, theta_hat] = least_sqr(u, num_samples, sysd);
writematrix(round(theta_hat,6));

%% ramp
num_samples = 500;

t = 0:num_samples-1;
unitstep = t>=0;
u = t.*unitstep;

[mse, sysd_predicted, theta_hat] = least_sqr(u, num_samples, sysd);
writematrix(round(theta_hat,6));

%% over parammeter

num_samples = 500;
input_variance = 1;
u = sqrt(input_variance) * randn(1, num_samples);
[mse, sysd_predicted, theta_hat] = least_sqr_over(u, num_samples, sysd);
writematrix(round(theta_hat,6));

%% under parammeter

num_samples = 500;
input_variance = 1;
u = sqrt(input_variance) * randn(1, num_samples);
[mse, sysd_predicted, theta_hat] = least_sqr_under(u, num_samples, sysd);
writematrix(round(theta_hat,6));

%% under parammeter LS


%% LS



