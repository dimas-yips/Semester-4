clear
clc
clf
%THT SISDIN 2023
%DIMAS YOGASTAMA PUTRA
%13621067

%% Inputs & Variables
NIM = input('Enter your NIM = ',"s"); %calculate values M, C, K, and w.
fprintf("\n");
M=1+str2double(NIM(6));
C0=3+str2double(NIM(4))+str2double(NIM(7));
C1=C0-0.5*C0;
C2=C0+0.5*C0;
K=6+str2double(NIM(5))+str2double(NIM(8));
w=0.5*(K/M)^0.5;
F0=1;

%% Analytical Solution
syms x1(t) x2(t)% define symbolic variables
eqns= [diff(x1,t)==x2, diff(x2,t)==(F0*sin(w*t)-K*x1-C0*x2)/M]; % Defines system of first-order ODE
conds=[x1(0)==0, x2(0)==0]; %initial conds
sols= dsolve (eqns, conds); %analytical solutions to the ODE system 
sol_x1 = sols.x1; %Analytical solution for 1st ODE
sol_x2 = sols.x2; %Analytical solution for 2nd ODE
fprintf('Analytical solution for 1st ODE: \n');
disp(sol_x1);
fprintf ('Analytical solution for 2nd ODE: \n');
disp(sol_x2);

%% Numerical Solution via Runge-Kutta
tspan= [0 25]; %specifies time interval 
X0 = [0 0];

% Define the values for C to loop over
Cvalues = [C0, C1, C2];
numCvalues = numel(Cvalues); %number of elements in the C_values array

% Plot the results for all values of C
figure(1)
subplot(2, 1, 1); % Create the first subplot for x(t)
hold on

for i = 1:numCvalues
    C = Cvalues(i);
    [tt, X] = ode45(@(tt,X)odefun(tt, X, M, C, K, F0, w), tspan, X0);
    
    plot(tt, X(:, 1), '-d', 'DisplayName', sprintf('solusi x(t) dengan C=%.2f', C));
end

title('Fungsi x(t)');
xlabel('t (s)');
ylabel('x (m)');
legend('show','FontSize', 7);

subplot(2, 1, 2); % Create the second subplot for dx(t)/dt
hold on

for i = 1:numCvalues
    C = Cvalues(i);
    [tt, X] = ode45(@(tt,X)odefun(tt, X, M, C, K, F0, w), tspan, X0);
    
    plot(tt, X(:, 2), '-^', 'DisplayName', sprintf('solusi dx(t)/dt dengan C=%.2f', C));
end

title('Fungsi dx(t)/dt');
xlabel('t (s)');
ylabel('dx/dt (m/s)');
legend('show','FontSize', 7)

% Plot the results for each value of C in separate subplots
figure(2)
for i = 1:numCvalues
    subplot(3, 1, i);
    C = Cvalues(i);
    [tt, X] = ode45(@(tt,X)odefun(tt, X, M, C, K, F0, w), tspan, X0);
    hold on
   
    plot(tt, X(:, 1), '-bd', 'DisplayName', sprintf('solusi x(t) dengan C=%.2f', C));
    plot(tt, X(:, 2), '-r^', 'DisplayName', sprintf('solusi dx(t)/dt dengan C=%.2f', C));
    ylim([-0.15, 0.15]);
    xlabel('t (s)');
    ylabel('x (m) & dx/dt (m/s)');
    legend('show','FontSize', 7);
    title(sprintf('Fungsi x(t) & dx(t)/dt dengan C=%.2f', C)); 
end
%% Function
function dXdt = odefun(tt, X, M, C, K, F, w)
    dXdt = zeros(2, 1);
    dXdt(1) = X(2);
    dXdt(2) = (F * sin(w * tt) - K * X(1) - C * X(2)) / M;
end