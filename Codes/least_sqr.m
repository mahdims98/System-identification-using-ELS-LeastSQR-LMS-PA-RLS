function [mse, sysd_predicted, theta_hat] = least_sqr(u, num_samples, sysd)
    rng(50)
    [num, den]=tfdata(sysd,'v');
    num_samples = num_samples;
    y_real = zeros(1, num_samples);
    phi = zeros(num_samples, 4+4);
    
    noise_variance = 0.01; % system noise variance
    noise = sqrt(noise_variance) * randn(1, num_samples);

    for i = 1:num_samples
        switch i
            case 1 % at time t=1
                phi_t = [zeros(1,4), zeros(1,4)];
            case 2
                phi_t = [-y_real(i-1),zeros(1,3), u(i-1),zeros(1,3)];
            case 3
                phi_t = [-y_real(i-1:-1:i-2),zeros(1,2), u(i-1:-1:i-2),zeros(1,2)];
            case 4
                phi_t = [-y_real(i-1:-1:i-3),zeros(1,1), u(i-1:-1:i-3),zeros(1,1)];
            otherwise
                phi_t = [-y_real(i-1:-1:i-4) u(i-1:-1:i-4)];
        end
        phi(i, :) = phi_t;
        y_real(i) = phi_t * [den(2:end), num(2:end)].' + noise(i);
    end
    
    %disp(det(phi.'*phi))
    Y_sample = y_real(1:num_samples).';
    theta_hat=(inv(phi.'*phi))*phi.'*Y_sample;
    y_predicted = zeros(1, num_samples);
    
%     for i = 1:num_samples
%         switch i
%             case 1 % at time t=1
%                 phi_t = [zeros(1,4), zeros(1,4)];
%             case 2
%                 phi_t = [-y_predicted(i-1),zeros(1,3), u(i-1),zeros(1,3)];
%             case 3
%                 phi_t = [-y_predicted(i-1:-1:i-2),zeros(1,2), u(i-1:-1:i-2),zeros(1,2)];
%             case 4
%                 phi_t = [-y_predicted(i-1:-1:i-3),zeros(1,1), u(i-1:-1:i-3),zeros(1,1)];
%             otherwise
%                 phi_t = [-y_predicted(i-1:-1:i-4) u(i-1:-1:i-4)];
%         end
%         y_predicted(i) = phi_t * theta_hat;
%     end

    y_predicted = phi * theta_hat;
    y_predicted = y_predicted.';


    mse = immse(y_real, y_predicted);

    % input
    figure();
    subplot(3,1,1);
    plot(u);
    title("Input signal");
    % y_real vs y_predicted
    
    subplot(3,1,2);
    plot(1:num_samples, y_real(1:num_samples), "--b", 'DisplayName','real output');
    xlabel("sample time");
    hold on;
    plot(1:num_samples, y_predicted(1:num_samples), 'DisplayName','predicted output');
    legend();
    title("Output")
  
    
    % error
    subplot(3,1,3);
    plot(1:num_samples, abs(y_predicted(1:num_samples)-y_real(1:num_samples)), 'DisplayName','error');
    xlabel("sample time")
    ylabel("abs. error")
    title("Abs. Error")
    
    % step response
    figure;
    sysd_predicted = tf([0, theta_hat(5:8).'], [1 , theta_hat(1:4).'], sysd.Ts);   
    bode(sysd_predicted);
    hold on; 
    bode(sysd);
    legend("Predicted system", "Real system");
end