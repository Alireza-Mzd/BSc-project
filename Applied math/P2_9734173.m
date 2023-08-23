%answer to the second question of Math Application Project
clear
clc

%defining constants
k1=0.0105;
k2=0.12;
t0=0;
tf=5; %seconds
h=0.001; %seconds
N=(tf-t0)/h;
t=t0:h:tf;

%defining Boundary Conditions
c0=0; %mol
c(1)=c0;
a0=110; %mol
a(1)=a0;
b0=220; %mol
b(1)=b0;

%defining ODEs
dadt=@(a,b,c) (k2*c-k1*a*b);
dbdt=@(a,b,c) (k2*c-k1*a*b);
dcdt=@(a,b,c) (3*k1*a*b-2*k2*c);

for i=1:N
    c(i+1)=c(i)+h*dcdt(a(i),b(i),c(i));
    b(i+1)=b(i)+h*dbdt(a(i),b(i),c(i));
    a(i+1)=a(i)+h*dadt(a(i),b(i),c(i));

end
fprintf('The final concentration of the substances a, b and c after 5 seconds are: %.4f, %.4f, %.4f moles respectively.',a(i+1),b(i+1),c(i+1))
plot(t,a,t,b,'r',t,c,'g')
legend('a','b','c','location','northeast')