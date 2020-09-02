%{ The Simple Gravity Pendulum}%
%Description: Plotting the angular position, angular velocity and phase
%space diagram for the simple pendulum using four methods: Forward Euler
%(linear and non-linear) and Trapezoidal (linear and non-linear)

clc, clear all, close all;
string = ['Select the Approximation Method:\n', ...
    '1. Forward Euler (linear ODE)\n',...
    '2. Forward Euler (non-linear ODE)\n',... 
    '3. Trapezoidal method (linear ODE)\n',...
    '4. Trapezoidal method (non-linear ODE)\n',...
    'Enter your selection\n'];


sel = input(string); % display the menu

%Initialize constant variables
g = 9.81; % the gravity of Earth
L = 1; % The length of the rod
ti= 0; % initial time
tf = 25; %final time
dt = 0.005; % time increment

nt = (tf-ti)/dt; 

t = linspace (ti, tf, nt);% create an array of time

theta = zeros(1, nt);% initialize the array of angular position 
omega = zeros(1, nt);% initialize the array of angular velocity

switch (sel)
    
    case 1
        %% Forward Euler (linear ODE)
        theta(1) = pi/10; % initial angular position for linear ODE
        
        for k = 1:nt-1
            omega(k+1) = omega(k) + dt*-g/L*theta(k);
            theta(k+1) = theta(k) + dt*omega(k);
        end 
      
        
    case 2
         %% Forward Euler (non-linear ODE)
         
        theta(1) = pi/2; % initial angular position for nonlinear ODE
         
        for k = 1:nt-1
            omega(k+1) = omega(k) + dt*-g/L*sin(theta(k));
            theta(k+1) = theta(k) + dt*omega(k);
        end
        
    
    case 3
        %% Trapezoidal method (linear ODE)
        
        theta(1) = pi/10; % initial angular position for linear ODE
        
        for k=1:nt-1
            A= [1 (dt/2)*g/L ; -dt/2 1];
            b = [omega(k)-dt/2*g/L*theta(k);dt/2*omega(k)+theta(k)];
            
            %Ax = b, so x[omega theta] = A\b
            x = A\b;
            omega(k+1) = x(1);
            theta(k+1) = x(2);
        end    
            
    
    case 4
        %% Trapezoidal method (non-linear ODE)  
        
        theta(1) = pi/2;  % initial angular position for nonlinear ODE
        err = 1e-15;
        
        d_theta = 10; %initialize the difference between new and old angular posititons
        d_omega = 10; %initialize the difference between new and old angular velocities
       
        
        for k = 1:nt-1
            %Using Forward Euler to get the initial values
            old_omega = omega(k) + dt*-g/L*sin(theta(k));
            old_theta = theta(k) + dt*omega(k);
            
            d_theta = 10; 
            d_omega = 10;
            while (d_theta > err && d_omega > err)
                %Using trapezoid method to calculate new values
                new_omega = omega(k) + dt/2*(-g/L*sin(theta(k)) - g/L*sin(old_theta));
                new_theta = theta(k) + (dt/2)*(omega(k) + old_omega);
                
                %Calculate the difference between new and old values
                d_omega = abs (new_omega - old_omega);
                d_theta = abs (new_theta - old_theta);
                
                old_omega = new_omega;
                old_theta = new_theta;
                
            
            end
            
            omega(k+1) = old_omega;
            theta(k+1) = old_theta;
        end
        
end    



%% Angular position graph

figure(1);
pos = plot(t, theta,'b');
xlabel ('time (s)');
ylabel ('theta (rad)');
title ('Angular position of the pendulum over time');
xlim ([0 25]);
ylim ([1.1*min(theta) 1.1*max(theta)]);
set(pos, 'LineWidth', 2);
set(gcf, 'Position', [100 40 1000 600]);
set(gca, 'LineWidth', 2, 'FontSize', 15);


%% Angular velocity graph

figure(2);
vel = plot (t, omega, 'r');
xlabel ('time (s)');
ylabel ('omega (rad/s)');
title ('Angular velocity of the pendulum over time');
xlim ([0 25]);
ylim ([1.1*min(omega) 1.1*max(omega)]);
set(vel, 'LineWidth', 2);
set(gcf, 'Position', [100 40 1000 600]);
set(gca, 'LineWidth', 2, 'FontSize', 15);



%% Phase space diagram
figure(3);
phase = plot (theta, omega);
xlabel ('theta (rad)');
ylabel ('omega (rad/s)');
title ('Phase space diagram of the pendulum');
xlim ([1.1*min(theta) 1.1*max(theta)]);
ylim ([1.1*min(omega) 1.1*max(omega)]);
set(phase, 'LineWidth', 2);
set(gcf, 'Position', [200 40 900 600]);
set(gca, 'LineWidth', 2, 'FontSize', 15);
