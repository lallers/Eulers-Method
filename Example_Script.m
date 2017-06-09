%%
clear all; clc                                     
tf = 70;                                          
h = 8.74e-4;
N = tf/h;
t = 0:h:tf;                                          

u(1,1) = 0;                      
u(2,1) = 1;  

for k = 1 : (length(t)-1) 
      u(:,k+1) = u(:,k) + h * euler3(t(k),u(:,k)); %Euler algorithm
end

u1ex = 1.582857756103056;
u2ex = -2.835763853877514;

for i = 1:length(t)
    Err(i) = sqrt(  (u(1,i) - u1ex)^2 + (u(2,i) - u2ex)^2   );
end

plot(t,Err)

u1 = u(1,end);
u2 = u(2,end);

%% Euler MEthod
clc;clear all;close all;
tf = 70;
h = [8.75e-4,8.75e-5,8.75e-6,8.75e-7];
uerr(1) = 1.582857756103056;
uerr(2) = -2.835763853877514;

for j = 1:length(h)

t = 0:h(j):tf;                                        
u = zeros(2,length(t));
u(1,1) = 0;                                          
u(2,1) = 1;

for i= 1 : (length(t)-1)                              
    u(:,i+1) = u(:,i) + h(j) * euler3(t(i),u(:,i));   
end

ErrE(j) = sqrt(  (u(1,end) - uerr(1))^2 + (u(2,end) - uerr(2))^2   );
hold on
plot(h(j),ErrE(j),'o')
legendInfo{j} = ['h = ' num2str(h(j))];
end

plot(h,ErrE,'b')

title('Error (Euler Method)');xlabel('Values of h');ylabel('Error(h)')
legend(legendInfo)
hold off
u1e = u(1,end);
u2e = u(2,end);

%%


tf = 70;
h = [.01,.001,.0001.00001.000001];

uerr(1) = 1.582857756103056;
uerr(2) = -2.835763853877514;

for j = 1:length(h)

t = 0:h(j):tf;                                        
u = zeros(2,length(t));
u(1,1) = 0;                                          
u(2,1) = 1;

for i= 1 : (length(t)-1)                              
    k_1 = euler3(t(i),u(:,i));
    k_2 = euler3(t(i)+0.5*h(j),u(:,i)+0.5*h(j)*k_1);
    k_3 = euler3((t(i)+0.5*h(j)),(u(:,i)+0.5*h(j)*k_2));
    k_4 = euler3((t(i)+h(j)),(u(:,i)+k_3*h(j)));

    u(:,i+1) = u(:,i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h(j);
    
end

ErrK(j) = sqrt(  (u(1,end) - uerr(1))^2 + (u(2,end) - uerr(2))^2   );
hold on
G = plot(h(j),ErrK(j),'o');
set(G, 'MarkerFaceColor', get(G, 'Color'));
legendInfo{j} = ['h = ' num2str(h(j))];

end

plot(h,ErrK,'b')
title('Error (4th - Runge Kutta)');xlabel('Values of h');ylabel('Error(h)')
legend(legendInfo)
hold off
u1k = u(1,end);
u2k = u(2,end);

%% Error Runtime
cutoff = 10e-5;
E = [ErrE(end) - cutoff;ErrK(end) - cutoff];

