%% MVCS LAB2
clc
close all
clear

%% Define matrices
alfa = 1.4;
beta = 1.8;
gamma = 1.05;
m = 13.5;
J = 0.37;
% State Matrix A (6x6)
A = [0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1;
     0 0 0 -alfa/m 0 0;
     0 0 0 0 -beta/m 0;
     0 0 0 0 0 -gamma/J];

% Thruster Allocation Matrix (3x4)
TAM = [0.7071 0.7071 -0.7071 -0.7071;    % Surge
       -0.7071 0.7071 -0.7071 0.7071;    % Sway
       -0.1888 0.1888 0.1888 -0.1888];   % Yaw

% Input Matrix B (6x4)
B_force_gain = [zeros(3,3);
                1/m 0 0;
                0 1/m 0;
                0 0 1/J];

K_gain = 30; % DC-gain from volts to force

B = B_force_gain * K_gain * TAM;

T = 1;
dt = 0.001;
xf = [-1 2 3 0 0 0]';

syms tau real
eA = expm(A*(T-tau)); % Symbolic matrix exponential
dG = eA * B * B' * eA'; % Gramian integrand
G = eval(int(dG,tau,0,T));

rankG = rank(G);
disp(['rank(G) = ' num2str(rankG)]);


if rankG == 6 % Check for full rank 
    disp('System is fully controllable.');
    %% Define time vector
    time = 0:dt:T;
    %time = time';
    LL = length(time);
    
    %% Compute system input
    iG = inv(G); % invert Gramian
    u(tau) = B'*eA'*iG*xf;
    u_time = u(time);
    uvalue = [eval(u_time{1}); eval(u_time{2}); eval(u_time{3}); eval(u_time{4})];    
    %% Initialize state
    x0 = [0 0 0 0 0 0]'; % we start from the origin
    x = zeros(6,LL); % this matrix will store the state evolution
    x(:,1) = x0;

    %% Simulate system (you can use Euler integration)
    for k = 1:LL-1
        u_k = uvalue(:,k); % Input vector at time k
        xdot = A * x(:,k) + B * u_k; % State derivative
        x(:,k+1) = x(:,k) + xdot * dt; % Euler integration step
    end
    
    %% Plot state evolution
    figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Set figure to full screen

    %% Plot inputs
    % Subplot 1: X Position
    subplot(2,3,1);
    plot(time, x(1,:),'linewidth',3);
    hold on;
    plot(time, ones(size(time)) * xf(1), 'linewidth',3);
    xlabel('Time[s]');
    ylabel('X[m]');
    legend('State','Reference');
    grid on;
    title('X');
    
    % Subplot 2: Y Position
    subplot(2,3,2);
    plot(time, x(2,:), 'linewidth',3);
    hold on;
    plot(time, ones(size(time)) * xf(2), 'linewidth',3);
    xlabel('Time[s]');
    ylabel('Y[m]');
    legend('State','Reference');
    grid on;
    title('Y');
    
    % Subplot 3: Theta
    subplot(2,3,3);
    plot(time, x(3,:), 'linewidth',3);
    hold on;
    plot(time, ones(size(time)) * xf(3), 'linewidth',3);
    xlabel('Time[s]');
    ylabel('Yaw[rad]');
    legend('State','Reference');
    grid on;
    title('Yaw');
    
    % Subplot 4: v_x
    subplot(2,3,4);
    plot(time, x(4,:), 'linewidth',3);
    hold on;
    plot(time, ones(size(time)) * xf(4), 'linewidth',3);
    xlabel('Time[s]');
    ylabel('v_x[m/s]');
    legend('State','Reference');
    grid on;
    title('v_x');
    
    % Subplot 5: v_y
    subplot(2,3,5);
    plot(time, x(5,:), 'linewidth',3);
    hold on;
    plot(time, ones(size(time)) * xf(5), 'linewidth',3);
    xlabel('Time[s]');
    ylabel('v_y[m/s]');
    legend('State','Reference');
    grid on;
    title('v_y');
    
    % Subplot 6: Yaw rate
    subplot(2,3,6);
    plot(time, x(6,:),'linewidth',3);
    hold on;
    plot(time, ones(size(time)) * xf(6), 'linewidth',3);
    xlabel('Time[s]');
    ylabel('Yaw rate[rad/s]');
    legend('State','Reference');
    grid on;
    title('Yaw rate');
    
    disp("final state - ref state = " + norm(x(:,end) - xf)); % should be small
else
   disp("System is not controllable.")
    [U,D,V] = svd(G); % SVD decomposition 
    disp('SVD analysis for 1.8:');
    disp('Singular values (D):');
    disp(diag(D));
    disp('Columns of U define the basis. Columns with non-zero singular values span the REACHABLE subspace.');
    disp('Columns with zero singular values span the UNREACHABLE subspace.');
end

%% Analysis 1.6: Thruster 2 Fails
disp('--- Analysis 1.6: Thruster 2 Fails ---');
TAM_mod1 = TAM;
TAM_mod1(:, 2) = 0; % Thruster 2 (second column) fails

B_mod1 = B_force_gain * K_gain * TAM_mod1;
% We need eA and tau from the base script
dG_mod1 = eA * B_mod1 * B_mod1' * eA';
G_mod1 = eval(int(dG_mod1, tau, 0, T));
rankG_mod1 = rank(G_mod1);
disp(['Rank with Thruster 2 failure = ' num2str(rankG_mod1)]);

if rankG_mod1 == 6
    disp('System IS STILL fully controllable.');
else
    disp('System is NOT controllable.');
end


%% Analysis 1.7: Thrusters 2 & 3 Fail
disp('--- Analysis 1.7: Thrusters 2 & 3 Fail ---');
TAM_mod2 = TAM;
TAM_mod2(:, [2, 3]) = 0; % Thrusters 2 and 3 fail

B_mod2 = B_force_gain * K_gain * TAM_mod2;
dG_mod2 = eA * B_mod2 * B_mod2' * eA';
G_mod2 = eval(int(dG_mod2, tau, 0, T));
rankG_mod2 = rank(G_mod2);
disp(['Rank with Thrusters 2 & 3 failure = ' num2str(rankG_mod2)]);

if rankG_mod2 < 6
    disp('System is NOT fully controllable.');
else
    disp('System IS STILL controllable.');
end