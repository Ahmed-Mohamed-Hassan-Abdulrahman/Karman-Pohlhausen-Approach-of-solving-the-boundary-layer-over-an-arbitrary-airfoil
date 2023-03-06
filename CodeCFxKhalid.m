clc
clear 
close all
%% Importing Airfoil Data
path="AirfoilData.xlsx";
path1="n0012-il.xlsx";
[X,U]=Import_Script(path);
[xfoil,yfoil]=Import_Script_airfoil(path1);
%% Givens
U_inf=50;
mu = 1.789*10^-5;
rho = 1.225;

%% Setting up The Paulhausen Method Functions
%
% $U',U'',\delta_1,\delta_2,\bar{\delta}_1,\bar{\delta}_2,K(\Lambda),f_1(\Lambda),f_2(\Lambda),Z$
U_d = zeros(1,length(U));
length_U = length(U);
%Velocity gradients
U_d=gradient(U,X);                                                         %U'
U_dd=gradient(U_d,X);                                                      %U''
%The Paulhausen Method Functions as Polynomials Coefficients
del1=[-120^-1 3/10];                                                       %delta1 [1st degree poly]
del2=[-1/9072 -1/945 37/315];                                              %delta2 [2nd degree poly]
k_p=conv(conv(del2,del2),[1 0]);                                           %K(Lambda) [5th degree poly]
f2_p=conv([1/6 2],del2);                                                   %f2(Lambda) [3rd degree]
F_del2=(conv((2*[0 0 f2_p]-4*k_p),del2)-[0 conv(2*k_p,del1)]);             %F(Lambda)*delta2 [7th degree poly]
lamda_s=roots(F_del2);                                                     %Solving for Lambda
a=real(lamda_s)<12 & real(lamda_s)>-12;                                    %Physical limits for Lambda

lamda_itr = zeros(1, length_U);
lamda_itr(1) = lamda_s(a);



delta1_ = zeros(1, length_U);
delta2_ = zeros(1, length_U);
delta1_(1) = 0.3-lamda_itr(1)/120;
delta2_(1) = 37/315-lamda_itr(1)/945-lamda_itr(1)^2/9072;

K = zeros(1, length_U);
K(1) = delta2_(1)^2*lamda_itr(1);

z = zeros(1, length_U);
z(1) = K(1)/U_d(1);


F = zeros(1, length_U);
f1 = zeros(1, length_U);
f2 = zeros(1, length_U);

f1(1) = delta1_(1)/delta2_(1);
f2(1) = (2+lamda_itr(1)/6)*delta2_(1);
F(1) = 2*f2(1)-4*K(1)-2*K(1)*f1(1);


F0_U0 = -0.0652*U_dd(1)/U_d(1)^2;
z(2) = z(1)+F0_U0*( X(2) - X(1) );

delta = zeros(1, length_U);
delta1 = zeros(1, length_U);
delta2 = zeros(1, length_U);

delta(1) = sqrt(lamda_itr(1)*mu/rho/U_d(1));
delta1(1) = delta(1)*delta1_(1);
delta2(1) = delta(1)*delta2_(1);


%% Populating The Paulhausen Method Functions

for i = 2:length(U_d)
    K(i) = z(i)*U_d(i);

     P1=[0 0 0 0 0 -K(i)]+k_p;
     Lamda_itr_s=roots(P1);
     Lamda_itr_s=Lamda_itr_s(imag(Lamda_itr_s)==0);
    b=real(Lamda_itr_s)<lamda_itr(1) & real(Lamda_itr_s)>-12;
    if b==0
        break
    end
    lamda_itr(i) = Lamda_itr_s(b);
    
    delta1_(i) = 0.3-lamda_itr(i)/120;
    delta2_(i) = 37/315-lamda_itr(i)/945-lamda_itr(i)^2/9072;
    
    delta(i) = sqrt(lamda_itr(i)*mu/rho/U_d(i));
    delta1(i) = delta(i)*delta1_(i);
    delta2(i) = delta(i)*delta2_(i);
    
    
    
    f1(i) = delta1_(i)/delta2_(i);
    f2(i) = (2+lamda_itr(i)/6)*delta2_(i);
    
    F(i) = 2*f2(i)-4*K(i)-2*K(i)*f1(i);
    
    z(i+1) = z(i)+F(i)/U(i)*(X(i+1)-X(i));
end
%% *Plots*
%% $\delta\:vs\:X_{Laminar}$
figure(1)
cla; hold on; grid on;  
set(gcf,'Color','White');                                              
set(gca,'FontSize',12);
plot(X, delta,'Color',[0.2 0.984 0.6],'LineWidth',2)
xlim([0 0.7])
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\delta$','interpreter','latex','FontSize',14);
title('$\delta$ $vs$ $X$','interpreter','latex','FontSize',14);
hold on
plot(xfoil(1:ceil(length(xfoil)/2)),yfoil(1:ceil(length(yfoil)/2))/12,'-.','Color',[0 0.6 0.96],'LineWidth',2)
legend('$\delta$ $vs$ $X_{Laminar}$','NACA0012','interpreter','latex','FontSize',8)

%% $\delta_1\:vs\:X_{Laminar}$
figure(2)
cla; hold on; grid on;  
set(gcf,'Color','White');                                               
set(gca,'FontSize',12);
plot(X, delta1,'Color',[0.2 0.984 0.6],'LineWidth',2)
xlim([0 0.7])
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\delta_{1}$','interpreter','latex','FontSize',14);
title('$\delta_{1}$ $vs$ $X$','interpreter','latex','FontSize',14);
hold on
plot(xfoil(1:ceil(length(xfoil)/2)),yfoil(1:ceil(length(yfoil)/2))/40,'-.','Color',[0 0.6 0.96],'LineWidth',2)
legend('$\delta_1$ $vs$ $X_{Laminar}$','NACA0012','interpreter','latex','FontSize',8)

%% $\delta_2\:vs\:X_{Laminar}$
figure(3)
cla; hold on; grid on;  
set(gcf,'Color','White');                                               
set(gca,'FontSize',12);
plot(X, delta2,'Color',[0.2 0.984 0.6],'LineWidth',2)
xlim([0 0.7])
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\delta_{2}$','interpreter','latex','FontSize',14);
title('$\delta_{2}$ $vs$ $X_{Laminar}$','interpreter','latex','FontSize',14);
hold on
plot(xfoil(1:ceil(length(xfoil)/2)),yfoil(1:ceil(length(yfoil)/2))/120,'-.','Color',[0 0.6 0.96],'LineWidth',2)
legend('$\delta_2$ $vs$ $X_{Laminar}$','NACA0012','interpreter','latex','FontSize',8)

%% $C_{fx}\:vs\:X_{Laminar}$
figure(4)
cla; hold on; grid on;  
set(gcf,'Color','White');                                               
set(gca,'FontSize',12);
cfx = ( 2.*f2)./(sqrt(rho*U_inf/mu).* U.* sqrt(z));
plot(X,cfx,'Color',[0.2 0.984 0.6],'LineWidth',2)
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$C{fx}$','interpreter','latex','FontSize',14);
title('$C_{fx}$ $vs$ $X_{Laminar}$','interpreter','latex','FontSize',14);
hold on
plot(xfoil(1:ceil(length(xfoil)/2)),yfoil(1:ceil(length(yfoil)/2))/12,'-.','Color',[0 0.6 0.96],'LineWidth',2)
legend('$C_{fx}$ $vs$ $X_{Laminar}$','NACA0012','interpreter','latex','FontSize',8)

%% $\tau_w\:vs\:X_{Laminar}$
figure(5)
cla; hold on; grid on;  
set(gcf,'Color','White');                                               
set(gca,'FontSize',12);
% tau_w = cfx * 0.5 * rho;
tau_w_notbar = 0.5*rho*cfx.*U.^2*U_inf^2;   
plot(X,tau_w_notbar,'Color',[0.2 0.984 0.6],'LineWidth',2)
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\tau_{w}$','interpreter','latex','FontSize',14);
title('$\tau_{w}$ $vs$ $X_{Laminar}$','interpreter','latex','FontSize',14);
hold on
plot(xfoil(1:ceil(length(xfoil)/2)),yfoil(1:ceil(length(yfoil)/2))*40,'-.','Color',[0 0.6 0.96],'LineWidth',2)
legend('$\tau_w$ $vs$ $X_{Laminar}$','NACA0012','interpreter','latex','FontSize',8)

%% Accounting for Transition Using Smith model to locate the point of transition.

H=f1;
log_Rex=log10(rho.*U.*X./mu);
F_H=-40.4557+64.8066*H-26.7538*H.^2+3.3819*H.^3;
d=find(abs(real(log_Rex-F_H))<0.01);
x_s=X(d);
f=f1(d);
delta2(d) = delta2(d-1);
delta(d) = 1.4*delta(d-1);

syms x
eqn= (-delta(d)/delta2(d)+x+3.3+1.5501*(x-0.6778)^-3.064==0);
H_(d) = double(vpasolve(eqn,x,1));
H_11 = double(vpasolve(1.1+0.86*(x-3.3)^-0.777==H_(d),x,1));


H_1(d)=H_11;    
for i = d+1:length(U_d)
    if H_(i-1) > 1.8 && H_(i-1) < 2.8
       break 
    end

    Re_delta2_(i-1) = U(i-1)*U_inf .* delta2(i-1) * (rho/mu);
    cfx(i-1) = 0.246*10^(-0.678*H_(i-1)) * (Re_delta2_(i-1))^-0.268;
    delta2_by_dx = cfx(i-1)/2 - (delta2(i-1) / U(i-1)) * U_d(i-1);
    delta2(i) = delta2(i-1) + delta2_by_dx * (X(i) - X(i-1));
    
    H_1(i) = (0.306e-1 * U(i) * ((H_1(i - 1) - 3) ^ (-0.6169e0)) + (U(i) * delta2(i) * H_1(i - 1) / (X(i) - X(i - 1)))) / ((U(i) * delta2(i) - U(i - 1) * delta2(i - 1)) / (X(i) - X(i - 1)) + U(i) * delta2(i) / (X(i) - X(i - 1)));
    
    if (H_1(i) <= 3.3)
        H_(i) = 3;
    elseif (H_1(i) > 5.3)
        H_(i) = 1.1+0.86*(H_1(i)-3.3)^-0.777;    
    else
        H_(i) = 0.6778 + 1.1536*(H_1(i) - 3.3 )^-0.326;
    end
end

%% $\delta_2\:vs\:X_{Transition}$
figure(6)

cla; hold on; grid on;  
set(gcf,'Color','White');                                               
set(gca,'FontSize',12);
plot(X, delta2,'Color',[0.2 0.984 0.6],'LineWidth',2)
xlim([0 1])
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\delta_{2}$','interpreter','latex','FontSize',14);
title('$\delta_{2}$ $vs$ $X$ $with$ $Transition$','interpreter','latex','FontSize',14);
hold on
plot(xfoil(1:ceil(length(xfoil)/2)),yfoil(1:ceil(length(yfoil)/2))/120,'-.','Color',[0 0.6 0.96],'LineWidth',2)
xline(x_s,'LineWidth',1.5)
legend('$\delta_2$ $vs$ $X_{Transition}$','NACA0012','Transition Line','interpreter','latex','FontSize',8)
%% $C_{fx}\:vs\:X_{Transition}$
figure(7)
cla; hold on; grid on;  
set(gcf,'Color','White');                                               
set(gca,'FontSize',12,'fontname','times');
plot(X,cfx,'Color',[0.2 0.984 0.6],'LineWidth',2)
hold on
plot(xfoil(1:ceil(length(xfoil)/2)),yfoil(1:ceil(length(yfoil)/2))/12,'-.','Color',[0 0.6 0.96],'LineWidth',2)
xline(x_s,'LineWidth',1.5)
legend('$C_{fx}$ $vs$ $X_{Transition}$','NACA0012','Transition Line','interpreter','latex','FontSize',8)
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$C{fx}$','interpreter','latex','FontSize',14);
title('$C_{fx}$ $vs$ $X$ $with$ $Transition$','interpreter','latex','FontSize',14);
%% $\tau_w\:vs\:X_{Transition}$ 
figure(8)
cla; hold on; grid on;                                                 
set(gcf,'Color','White');                                             
set(gca,'FontSize',12,'fontname','times');
tau_w = cfx * 0.5 * rho.*U.^2*U_inf^2;
plot(X,tau_w,'Color',[0.2 0.984 0.6],'LineWidth',2)
hold on
plot(xfoil(1:ceil(length(xfoil)/2)),yfoil(1:ceil(length(yfoil)/2))*40,'-.','Color',[0 0.6 0.96],'LineWidth',2)
xline(x_s,'LineWidth',1.5)
legend('$\tau_{w}$ $vs$ $X_{Transition}$','NACA0012','Transition Line','interpreter','latex','FontSize',8)
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\tau_{w}$','interpreter','latex','FontSize',14);
title('$\tau_{w}$ $vs$ $X$ $with$ $Transition$','interpreter','latex','FontSize',14);
%%
X_req=linspace(0.1,1,10);
index=zeros(1,length(X_req));
U_req=zeros(1,length(X_req));
U_d_req=zeros(1,length(X_req));
Z_req=zeros(1,length(X_req));
K_req=zeros(1,length(X_req));
lamda_itr_req=zeros(1,length(X_req));
f1_req=zeros(1,length(X_req));
f2_req=zeros(1,length(X_req));
F_req=zeros(1,length(X_req));
delta_req=zeros(1,length(X_req));
delta1_req=zeros(1,length(X_req));
delta2_req=zeros(1,length(X_req));
cfx_req=zeros(1,length(X_req));
tau_w_notbar_req=zeros(1,length(X_req));
for i = 1:1:length(X_req)  
index(i)=find(X>=X_req(i),1,'first' );
U_req(i)=U(index(i))*U_inf;
U_d_req(i)=U_d(index(i))*U_inf;
Z_req(i)=z(index(i));
K_req(i)=K(index(i));
lamda_itr_req(i)=lamda_itr(index(i));
f1_req(i)=f1(index(i));
f2_req(i)=f2(index(i));
F_req(i)=F(index(i));
delta_req(i)=delta(index(i));
delta1_req(i)=delta1(index(i));
delta2_req(i)=delta2(index(i));
cfx_req(i)=cfx(index(i));
tau_w_notbar_req(i)=tau_w_notbar(index(i));
end
%% Getting Required Data to Fill the Table
fprintf('======= Sussy results =======\n');
Data = [X_req',U_req',U_d_req',Z_req',K_req',lamda_itr_req',f1_req',f2_req',delta_req',delta1_req',delta2_req',cfx_req',tau_w_notbar_req'];
VarNames = {'X','U','U_dash','Z','K','Lambda','f1','f2','delta','delta 1','delta 2','Cfx','tau_w'};
T = table(Data(:,1),Data(:,2),Data(:,3),Data(:,4),Data(:,5),Data(:,6),Data(:,7),Data(:,8),Data(:,9),Data(:,10),Data(:,11),Data(:,12),Data(:,13),'VariableNames',VarNames)
   warning('off','MATLAB:xlswrite:AddSheet'); %optional
   writematrix(Data,'test.xlsx','Sheet',1);