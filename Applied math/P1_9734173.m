%Answer to the first question of Math Application Project
clear
clc
%defining parameters
C=1.275;
Do=2.75; %cm
A=3.14/4*(Do*10^(-2))^2; %m2
g=9.81; %m/s2
R=3; %m
%ODE
dHdt=@(H) (-(C*A*sqrt(2*g))/(pi*(2*R*sqrt(H)-H^(1.5))));
%Boundary Condition
H(1)=3.75;
Hf=0; %final height, m
t0=0; %s
h=1; 
i=1;
 while H(i)>0

     H(i+1)=H(i)+dHdt(H(i))*h;
     i=i+1;
     
 end
height=H(i-1);
format long
j=((i-1)*h); %time required, seconds
fprintf('-the required time to drain the tank is about %.5f seconds  or %.5f minutes or %.5f hours. ' , j, j/60, j/3600)

%plotting part :
time=(0:length(H)-1)/60;
plot(time,H,'g')
title('Height of Fluid Vs. Time')
xlabel('Time (Minutes)')
ylabel('Height (m)')
axis([0 300 0 4])
%---------------------------------------------------------------------------
%as an EXTRA, by calculating the real time using analytical method we have:
syms H
I=(1/(-(C*A*sqrt(2*g))/(pi*(2*R*sqrt(H)-H^(1.5)))));
a=int(I,H,3.75,0);
a=a/60;
b=abs((double(a)-j/60))/a*100;
%---------------------------------------------------------------------------
fprintf('\n-the real required time is %.5f minutes and therefore the Error is %.3f percent.',a,b)
fprintf('\n-the final height is %s meters which is not exactly zero,\n  but since it is obtained by a numerical method it can be considered as zero.',height)
