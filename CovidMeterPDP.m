clear;
clc;

%% Research Code by Agus Hasan

load DATA.mat;  % date | month | susceptible cases (S) | probable cases (P) | active cases (I) | recovered cases (R) | death cases (D)
load RtW.mat;   % Rt without PDP

Tinf     = 12;      % infectious time
Std_Tinf = 3;       % standard deviation of infectious time
sigma    = 1.96;    % 95% confidence interval
N        = sum(DATA(1,3:end));
tf       = length(DATA(:,1));
td       = datetime(2020,DATA(1,2),DATA(1,1)-1) + caldays(1:tf);
CFR      = DATA(end,end)/(sum(DATA(end,5:end)));
dt       = 0.01;
t        = dt:dt:tf;
x2       = [td, fliplr(td)];

%% Measurement/data matrix
C = [1 0 0 0 0 0;
     0 1 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 1 0];
 
%% Noise
QF = diag([10 10 10 5 5 0.2]);   % process and measurement covariance matrices
RF = diag([100 10 5 5 1]);       % are considered as tuning parameters

%% Numerical Simulation
for m = 1:3
% Infection time
Ti     = Tinf-Std_Tinf*sigma+(m-1)*Std_Tinf*sigma;    % infection time with standard dev. 1 day
gamma  = (1-CFR)/Ti;    % recovery date
delta  = CFR/Ti;        % death rate
omega  = delta;         % death rate
eps    = gamma+delta;   % negative testing rate
kappa  = 0.08;          % positive testing rate

%% Initialization
xhat     = [N-2; 1; 1; 0; 0; 1];        % initial condition
xhatBeta = 0;                           
Pplus    = 1000000*eye(length(xhat));
% for plotting
xhatArray    = [];

%% Extended Kalman Filter
for i=1:((tf-1)/dt)
     xhatArray    = [xhatArray xhat];    
     % assimilating confirmaed cases
     y = [interp1(0:1:tf-1,DATA(:,3),t);    % S
         interp1(0:1:tf-1,DATA(:,4),t);     % P
         interp1(0:1:tf-1,DATA(:,5),t);     % I
         interp1(0:1:tf-1,DATA(:,6),t);     % R
         interp1(0:1:tf-1,DATA(:,7),t)];    % D
     % simulating the model
     xhat(1) = xhat(1)-xhat(6)*xhat(1)*(xhat(2)+xhat(3))*dt/N+eps*dt*xhat(2);
     xhat(2) = xhat(2)+xhat(6)*xhat(1)*xhat(2)*dt/N-(eps+kappa+omega)*xhat(2)*dt;
     xhat(3) = xhat(3)+xhat(6)*xhat(1)*xhat(3)*dt/N+kappa*xhat(2)*dt-(gamma+delta)*xhat(3)*dt;
     xhat(4) = xhat(4)+gamma*xhat(3)*dt;
     xhat(5) = xhat(5)+delta*xhat(3)*dt+omega*xhat(2)*dt;
     xhat(6) = xhat(6);
    % calculating the Jacobian matrix
    FX    = [1-xhat(6)*(xhat(2)+xhat(3))*dt/N -xhat(6)*xhat(1)*dt/N+eps*dt -xhat(6)*xhat(1)*dt/N 0 0 -xhat(1)*(xhat(2)+xhat(3))*dt/N;
             xhat(6)*xhat(2)*dt/N 1+xhat(6)*xhat(1)*dt/N-(eps+kappa+omega)*dt 0 0 0 xhat(1)*xhat(2)*dt/N;
             xhat(6)*xhat(3)*dt/N kappa*dt 1+xhat(6)*xhat(1)*dt/N-(gamma+delta)*dt 0 0 xhat(1)*xhat(3)*dt/N;
             0 0 gamma*dt 1 0 0;
             0 omega*dt delta*dt 0 1 0;
             0 0 0 0 0 1];
    Pmin  = FX*Pplus*FX'+QF;
    % update 
    KF    = Pmin*C'*inv(C*Pmin*C'+RF);  % Kalman gain
    xhat  = xhat + KF*(y(:,i)-C*xhat);
    Pplus = (eye(6)-KF*C)*Pmin;
    xhat(6) = max(0,(xhat(1)/N)*xhat(6));
end

windowSize = 300;
b          = (1/windowSize)*ones(1,windowSize);
a          = 1;
xhatArray(5,:) = filter(b,a,xhatArray(5,:));

xhatSArray  = [];
xhatS       = xhatArray(1,tf);
xhatPArray  = [];
xhatP       = xhatArray(2,tf);
xhatIArray  = [];
xhatI       = xhatArray(3,tf);
xhatRArray  = [];
xhatR       = xhatArray(4,tf);
xhatDArray  = [];
xhatD       = xhatArray(5,tf);
xhatBArray  = [];
xhatB       = xhatArray(6,tf);
for i=1:tf-1
    xhatSArray  = [xhatSArray xhatS];
    xhatS       = xhatArray(1,(1/dt)*i);
    xhatPArray  = [xhatPArray xhatP];
    xhatP       = xhatArray(2,(1/dt)*i);
    xhatIArray  = [xhatIArray xhatI];
    xhatI       = xhatArray(3,(1/dt)*i);
    xhatRArray  = [xhatRArray xhatR];
    xhatR       = xhatArray(4,(1/dt)*i);
    xhatDArray  = [xhatDArray xhatD];
    xhatD       = xhatArray(5,(1/dt)*i);
    xhatBArray  = [xhatBArray xhatB];
    xhatB       = xhatArray(6,(1/dt)*i);
end

xhatSArray  = [xhatSArray xhatS];
xhatPArray  = [xhatPArray xhatP];
xhatIArray  = [xhatIArray xhatI];
xhatRArray  = [xhatRArray xhatR];
xhatDArray  = [xhatDArray xhatD];
xhatBArray  = [xhatBArray xhatB];

%% Calculating the effective reproduction number
R0Array = [];
R0 = 0;
for j = 1:tf
    R0Array = [R0Array R0];
    M = [xhatBArray(j)/(eps+kappa+omega) 0;
         xhatBArray(j)*kappa/((gamma+delta)*(eps+kappa+omega)) xhatBArray(j)/(gamma+delta)];
    R0 = max(xhatBArray(j)/(eps+kappa+omega),xhatBArray(j)/(gamma+delta));
end

Z(m,:) = R0Array;

end

for l = 1:tf
    curve2(l)      = max(Z(:,l));
    Rt(l)          = mean(Z(:,l));
    curve1(l)      = min(Z(:,l));
end

curve1 = smooth(curve1);
curve2 = smooth(curve2);

figure(1)
subplot(2,2,1)
plot(td,DATA(:,4),'ob','LineWidth',6)
hold on
plot(td,xhatPArray,'-r','LineWidth',4);
xlim([min(td) max(td)])
ylabel('Probable Cases')
set(gca,'color','none','FontSize',36)
legend('Actual','Estimated')
grid on;
grid minor;
subplot(2,2,2)
plot(td,DATA(:,5),'ob','LineWidth',6);
hold on
plot(td,xhatIArray,'-r','LineWidth',4)
xlim([min(td) max(td)])
ylabel('Active Cases')
set(gca,'color','none','FontSize',36)
legend('Actual','Estimated')
grid on;
grid minor;
subplot(2,2,3)
plot(td,DATA(:,6),'ob','LineWidth',6);
hold on
plot(td,xhatRArray,'-r','LineWidth',4)
xlim([min(td) max(td)])
ylabel('Recovered Cases')
set(gca,'color','none','FontSize',36)
legend('Actual','Estimated')
grid on;
grid minor;
subplot(2,2,4)
plot(td,DATA(:,7),'ob','LineWidth',6);
hold on
plot(td,xhatDArray,'-r','LineWidth',4)
xlim([min(td) max(td)])
ylabel('Death Cases')
set(gca,'color','none','FontSize',36)
legend('Actual','Estimated')
alpha(0.3)
grid on;
grid minor;

figure(2)
yyaxis left
bar(td,DATA(:,5),'g')
hold on
bar(td,DATA(:,4),'y')
ylabel('Cases','color','k')
alpha(0.1)
yyaxis right
inBetween = [curve1', fliplr(curve2')];
fill(x2, inBetween, 'b');
hold on;
plot(td,smooth(Rt),'-b','LineWidth',4);
hold on
inBetween = [RtW(1,:), fliplr(RtW(3,:))];
fill(x2, inBetween, 'r');
hold on;
plot(td,RtW(2,:),'-r','LineWidth',4)
hold on;
plot(td,ones(1,tf),':c','LineWidth',4)
set(gca,'color','none','FontSize',36)
ylabel('Rt','color','k')
xlim([min(td) max(td)])
grid on;
grid minor;
alpha(0.3)
legend('Active Cases','Probable Cases','95% CI with Probable Cases','Rt with Probable Cases','95% CI without Probable Cases','Rt without Probable Cases','Rt =1','FontSize',24)
ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';

figure(3)
STC = [DATA(:,4) DATA(:,5) DATA(:,6) DATA(:,7)];
bar(td,STC,'stacked')
set(gca,'color','none','FontSize',36)
grid on;
grid minor;
alpha(0.8)
legend('Probable Cases','Active Cases','Recovered Cases','Death Cases')