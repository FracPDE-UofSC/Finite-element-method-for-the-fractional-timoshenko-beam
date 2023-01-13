%run timo
clear
clc
E_alp=1.99948e11;
height=0.025;
Area=0.1*height; I=0.1*height^3/12;
rho=8193.2518;
ks=(5/6)^2;
G_alp=7.99792e10;
L=1;
T=0.1;
numberElements=2048;
f=0;
condition='healthy';
m1=1;
bd='CC';
alpha=0.5;
omega=timo_fre(1,height,E_alp,G_alp,rho,ks,L);
[w2,w3,w_euler]=getfre(E_alp,G_alp,Area,I,rho,ks,L,1);
M0=64;
tau=2*pi/omega/M0;
M=ceil(T/tau);
h=1/numberElements;
x_points=0:32*h:L;
time=0:32*tau:M*tau;
[w1_alp,phi_alp]=fractimoL1(T,M0,numberElements,alpha,height,omega);
