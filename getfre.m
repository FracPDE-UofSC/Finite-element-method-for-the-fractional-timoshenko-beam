function [w2,w3,w1]=getfre(E,G,Area,I,rho,ks,l,n)
%% get natural frequency for SS Timoshenko and Euler
alp=sqrt(E*I/rho/Area);
r=sqrt(I/Area);
a1=rho*r^2/ks/G;
a2=1+n^2*pi^2*r^2/l^2+n^2*pi^2*r^2/l^2*E/ks/G;
a3=alp^2*n^4*pi^4/l^4;
%     syms a1 a2 a3 x
%     equ=a1*x.^4-a2*x.^2+a3==0;
%     solve(equ,x)


%Euler SS
w1=sqrt(alp^2*n^4*4.73004074^4/l^4);
%Timoshenko SS
w2= ((a2 - (a2^2 - 4*a1*a3)^(1/2))/(2*a1))^(1/2);
w3=((a2 + (a2^2 - 4*a1*a3)^(1/2))/(2*a1))^(1/2);

end