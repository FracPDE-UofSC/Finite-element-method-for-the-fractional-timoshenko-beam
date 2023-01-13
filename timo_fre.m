function omega_timo=timo_fre(n,h,E,G,rho,k,L)
% n=1;
% h=0.1;
height=h;
% E=1.99948e11;
Area=0.1*height; I=0.1*height^3/12;
% L=1;
% rho=8193.2518;
% G=7.99792e10;
b=@(p)sqrt(rho*Area*L^4*p^2/E/I);
% k=(5/6)^2;% G=8.21e8;
r=sqrt(I/Area/L^2);
s=sqrt(E*I/k/Area/G/L^2);
alpha=@(p)sqrt(-(r^2+s^2)+sqrt((r^2-s^2)^2+4/b(p)^2))/sqrt(2);
beta=@(p)sqrt((r^2+s^2)+sqrt((r^2-s^2)^2+4/b(p)^2))/sqrt(2);
% f=@(p)sin(b(p)*beta(p));%ss
f=@(p)2-2*cosh(b(p)*alpha(p))*cos(b(p)*beta(p))+b(p)/sqrt(1-b(p)^2*r^2*s^2)...
    *(b(p)^2*s^2*(r^2-s^2).^2+(3*s^2-r^2))*sinh(b(p)*alpha(p))*sin(b(p)*beta(p));
	omega_timo = fsolve(f,700); %300
	if(h > 0.025)
omega_timo = fsolve(f,3000);
end

%% SS
% n=1;
% alpha2=E*I/rho/Area;
% r2=I/Area;
% a1=rho*r2/k/G;
% a2=(1+n^2*pi^2*r2/L^2+n^2*pi^2*r2/L^2*E/k/G);
% a3=alpha2*n^4*pi^4/L^4;
% f2=@(w) a1*w.^4-a2*w^2+a3;
% omega_timo = fsolve(f2,700);
% omega_euler=sqrt(alpha2*n^4*pi^4/L^4);
% end