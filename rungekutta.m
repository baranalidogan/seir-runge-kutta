%First, run modelcovidSEIR.m, then the rungekutta.m
% Initial Conditions and Simulation Time
y0 = [7610026000, 80000, 24545, 907]; % y0 = [v0, w0]
t = linspace(0,60,601)';
% Runge-Kutta 4th-Order Algorithm
y_Kutta = zeros(length(t), 4);
y_Kutta(1, :) = y0;
h = t(2)-t(1); % Constant time step
for i = 2:length(t)
k1 = modelcovidSEIR(t(i-1), y_Kutta(i-1, :));
k2 = modelcovidSEIR(t(i-1)+h/2, y_Kutta(i-1, :)+k1*h/2);
k3 = modelcovidSEIR(t(i-1)+h/2, y_Kutta(i-1, :)+k2*h/2);
k4 = modelcovidSEIR(t(i-1)+h, y_Kutta(i-1, :)+k3*h);
y_Kutta(i, :) = y_Kutta(i-1, :)+(k1/6+k2/3+k3/3+k4/6)*h;
end
figure; subplot(4,1,1); hold on;
plot(t,y_Kutta(:,1),'blue-s','MarkerIndices',1:50:length(t))
%plot(t, y_Kutta(:,1),'o',t, uapprox, 'red');
grid on; hold off
%title('RK4', 'Interpreter', 'Latex')
subplot(4,1,2); hold on;
plot(t,y_Kutta(:,2),'blue-s','MarkerIndices',1:50:length(t))
%plot(t, y_Kutta(:,2),'o',t, wapprox, 'red');
grid on; hold off
subplot(4,1,3); hold on;
plot(t,y_Kutta(:,3),'blue-s','MarkerIndices',1:50:length(t))
grid on; hold off
subplot(4,1,4); hold on;
plot(t,y_Kutta(:,4),'blue-s','MarkerIndices',1:50:length(t))
grid on; hold off