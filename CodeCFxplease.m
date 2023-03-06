clc
clear
close all

%% Importing the airfoil characteristics
path="AirfoilData.xlsx";
[X,U]=Import_Script(path);
U_d = zeros(1,length(U));                                                  %U' intialization
length_U = length(U);

U_d=gradient(U,X);                                                         %U' calculation
U_dd=gradient(U_d,X);                                                      %U'' calculation
Uinf=50;

del1=[-120^-1 3/10];                                                       %delta1 polynomial representation in pohlhausen parameter *1st degree*
del2=[-1/9072 -1/945 37/315];                                              %delta2 polynomial representation in pohlhausen parameter *2nd degree*
k_p=conv(conv(del2,del2),[1 0]);                                           %K=delta2^2*(Pohlhausen parameter) polynomial  *5th order*
f2_p=conv([1/6 2],del2);                                                   %f2 in polynomial form = delta2*(2+(pohl parameter)/6) *3rd degree*
F_del2=(conv((2*[0 0 f2_p]-4*k_p),del2)-[0 conv(2*k_p,del1)]);             %F*delta2=7agat (polynomial) %zeros are to adjust to 5th degree, first term will be a 7th degree,
lamda_s=roots(F_del2);                                                     %roots of the polynomial F*delta2
a=real(lamda_s)<12 & real(lamda_s)>-12;                                    %checking which root is within the range -12<lambda<12

lamda_itr = zeros(1, length_U);                                            %Lambda for each position/ velocity on the airfoil
lamda_itr(1) = lamda_s(a);                                                 %

mu = 1.789*10^-5;
rho = 1.225;

delta1_ = zeros(1, length_U);
delta2_ = zeros(1, length_U);
delta1_(1) = 0.3-lamda_itr(1)/120;                                         %Delta 1 (non_dimensionalized)
delta2_(1) = 37/315-lamda_itr(1)/945-lamda_itr(1)^2/9072;                  %Delta 2 (non_dimensionalized)

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
% 



for i = 2:length(U_d)
    K(i) = z(i)*U_d(i);

     P1=[0 0 0 0 0 -K(i)]+k_p;
     Lamda_itr_s=roots(P1);
    j=abs(imag(Lamda_itr_s))>0;
    Lamda_itr_s(j)=[];
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

% Delta vs X
figure(1)
plot(X, delta)
xlim([0 0.7])
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\delta$','interpreter','latex','FontSize',14);
title('$\delta$ $vs$ $X$','interpreter','latex','FontSize',14);


% Delta1 vs X
figure(2)
plot(X, delta1)
xlim([0 0.7])
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\delta_{1}$','interpreter','latex','FontSize',14);
title('$\delta_{1}$ $vs$ $X$','interpreter','latex','FontSize',14);


% Delta2 vs X
figure(3)
plot(X, delta2)
xlim([0 0.7])
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\delta_{2}$','interpreter','latex','FontSize',14);
title('$\delta_{2}$ $vs$ $X$ $Laminar$','interpreter','latex','FontSize',14);


% Cfx vs X
figure(4)
cfx = ( 2.*f2 )./ (sqrt(rho*Uinf./mu).*U.* sqrt(z));
cfx_bar=2/(rho*Uinf./mu)*sqrt(rho*Uinf/mu)*f2.*U./sqrt(z);
plot(X,cfx)
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$C{fx}$','interpreter','latex','FontSize',14);
title('$C_{fx}$ $vs$ $X$ $Laminar$','interpreter','latex','FontSize',14);

% tau_w vs x
tau_w_notbar = 0.5*rho*cfx.*U.^2*Uinf^2;
figure(5)
plot(X,tau_w_notbar)
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\tau_{w}$','interpreter','latex','FontSize',14);
title('$\tau_{w}$ $vs$ $X$ $Laminar$','interpreter','latex','FontSize',14);


%% Accounting for Transition

H=f1;
log_Rex=log10(rho.*U.*X./mu);
F_H=-40.4557+64.8066*H-26.7538*H.^2+3.3819*H.^3;
d=find(abs(log_Rex-F_H)<0.01);
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

    Re_delta2_(i-1) = Uinf*U(i-1) .* delta2(i-1) * (rho/mu);
    cfx(i-1) = 0.246*10^(-0.678*H_(i-1)) * (Re_delta2_(i-1))^-0.268;
    delta2_by_dx = cfx(i-1)/2 - ((delta2(i-1) / U(i-1)) * U_d(i-1))*(H_(i-1)+2);
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


% Delta2 vs X
figure(6)
plot(X, delta2)
xlim([0 1])
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\delta_{2}$','interpreter','latex','FontSize',14);
title('$\delta_{2}$ $vs$ $X$ $with Transition$','interpreter','latex','FontSize',14);


% Cfx vs X
figure(7)
plot(X,cfx)
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$C{fx}$','interpreter','latex','FontSize',14);
title('$C_{fx}$ $vs$ $X$ $with Transition$','interpreter','latex','FontSize',14);

% tau_w vs x 
tau_w = cfx * 0.5 * rho.*U.^2*Uinf^2;
figure(8)
plot(X,tau_w)
xlabel('$X$','interpreter','latex','FontSize',14);
ylabel('$\tau_{w}$','interpreter','latex','FontSize',14);
title('$\tau_{w}$ $vs$ $X$ $Transition$','interpreter','latex','FontSize',14);

