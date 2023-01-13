function [w1,phi1]=fractimoL1(T,M,numberElements,alpha,h,omega)
% M=100;
% numberElements=48;
condition='healthy';
m1=1;
bd='CC';
height=h;
% E_alp=3.4e9;
% rho=1.22e3;
% ks=0.35;
% G_alp=1.2e9;
%Area=0.1*height; I=0.1*height^3/12;

E_alp=1;
rho=1;
ks=1;
G_alp=1;
Area=1; I=1;



% E_alp=1.99948e11;
% rho=8193.2518;
% ks=(5/6)^2;
% G_alp=7.99792e10;

amplitude=1;
points(1).x=0;
points(1).y=0;
points(2).x=1;
points(2).y=0;
L=points(2).x-points(1).x;
n0=1; %1st natural frequency
r=1;
time=zeros(M+1,1);
for i=1:M+1
    time(i)=T*((i-1)/M)^r;
end
tau=diff(time);
Iwt=@(t)gamma(2.8)/gamma(3.8-alpha)*t.^(2.8-alpha);
wtt=@(t)1.8*t.^(0.8);
% wx=@(x)x.^2.*(1-x).^2;
% wxx=@(x)4.*x.*(2*x - 2) + 2*(x - 1).^2 + 2*x.^2;
% wxxx=@(x)2-12*x+12*x.^2;
wx=@(x)sin(pi*x);
wxx=@(x)pi*cos(pi*x);
wxxx=@(x)-pi^2*sin(pi*x);
%%  BC
%%pined-pinned SSBC
if bd=='SS'
fixedNodeU =[1 numberElements+1]';
%  fixedNodeV=[1 numberElements+1]';
fixedNodeV=[]';
end
if bd=='CC'
  fixedNodeU =[1 numberElements+1]';
  fixedNodeV=[1 numberElements+1]';
end
if bd=='CS'
    fixedNodeU =[1 numberElements+1]';
    fixedNodeV=[1 ]';
end

if bd=='CT'
    fixedNodeU =[1 ]';
    fixedNodeV=[1 ]';
end
%%fixed-free
% fixedNodeU =[1]';
% fixedNodeV =[2]';

numberNodes=numberElements+1;
h=L/numberElements;
x_nodes=0:h:L;
GDof=numberNodes;
activeDof_w=setdiff([1:GDof]',[fixedNodeU]);
activeDof_phi=setdiff([1:GDof]',[fixedNodeV]);
e=ones(numberElements+1,1);
if condition=='healthy'
    Stiff=spdiags([-1/h*e 2/h*e -1/h*e], -1:1,numberElements+1, numberElements+1 );
end
if condition=='damaged'
    Stiff=damaged_stiff(points(1).x,points(2).x,numberElements);
end
Stiff(1,1)=1/h;
Stiff(end,end)=1/h;
Mass=spdiags([h*e/6 2*h*e/3 h*e/6], -1:1,numberElements+1,numberElements+1  );
Mass(1,1)=h/3;
Mass(end,end)=h/3;
A=spdiags([-0.5*e 0*e 0.5*e], -1:1,numberElements+1, numberElements+1 );
A(1,1)=-1/2;
A(end,end)=1/2;
b=zeros(M+1,1);
w=zeros(length(activeDof_w),M+1);
phi=zeros(length(activeDof_phi),M+1);
% omega=(beta^2*sqrt(E_alp*I/rho/Area));
fq=@(t)amplitude*cos(omega*t);
% intq=@(t)amplitude*sin(omega*t)/omega;
intq=@(t)0.01*t;
u0=@(x)1;
% rhsb=h*ones(numberElements+1,1);%distributed load
% rhsb(1)=h/2;
% rhsb(end)=h/2;
rhsb=rhs_liner(u0,numberElements,L,activeDof_w); %Xn

% rhsb=zeros(numberElements+1,1); %point load
% rhsb(40,1)=2;
for n=1:M %2:M
     tn=time(n+1);
    for k=1:n-1%n-1
%             b(k)=((tn-time(k))^(1-alpha)-(tn-time(k+1))^(1-alpha))/gamma(2-alpha);
        b(k)=((tn-time(k))^(2-alpha)-2*(tn-time(k+1))^(2-alpha)+(tn-time(k+2))^(2-alpha) )/tau(n)/gamma(3-alpha);
    end
    b(n)=tau(n)^(1-alpha)/gamma(3-alpha);
    %% integral
    A1=rho*Area*Mass(activeDof_w,activeDof_w)...
        +ks*tau(n)*G_alp*Area*b(n)*Stiff(activeDof_w,activeDof_w);
    b1=rhsb*tau(n)*(intq(tn)-intq(0));
    b1=b1+rho*Area*Mass(activeDof_w,activeDof_w)*w(:,n);
    for k=1:n-1
        b1=b1...
            -tau(n)*ks*G_alp*Area*b(k)*Stiff(activeDof_w,activeDof_w)*w(:,k+1);
    end
    
    B1= tau(n)*ks*G_alp*Area*b(n)*A(activeDof_phi,activeDof_phi);
    for k=1:n-1
        c=tau(n)*ks*G_alp*Area*b(k)*A(activeDof_phi,activeDof_phi)*phi(:,k+1);
        if bd=='SS'
        b1=b1-c(activeDof_w);
        else
           b1=b1-c;
        end
    end
    B2=rho*I*Mass(activeDof_phi,activeDof_phi)...
        +  tau(n)*E_alp*I*b(n)*Stiff(activeDof_phi,activeDof_phi)...
        + tau(n)*b(n)*ks*G_alp*Area*Mass(activeDof_phi,activeDof_phi);
    b2=rho*I*Mass(activeDof_phi,activeDof_phi)*phi(:,n);
    for k=1:n-1
        b2=b2...
            - tau(n)*E_alp*I*b(k)*Stiff(activeDof_phi,activeDof_phi)*phi(:,k+1)...
            -tau(n)*b(k)*ks*G_alp*Area*Mass(activeDof_phi,activeDof_phi)*phi(:,k+1);
    end
    B3= - tau(n)*b(n)*ks*G_alp*Area*A(activeDof_phi,activeDof_phi);
    for k=1:n-1
        if bd=='SS'
        d=[0;w(:,k+1);0];
        else
        d=w(:,k+1);
        end   
        b2=b2+tau(n)*b(k)*ks*G_alp*Area*A(activeDof_phi,activeDof_phi)*d;
    end
    AA=[A1, B1; B3, B2]; %CC
    bb=[b1;b2];
    x_sol=AA\bb;
   w(:,n+1)=x_sol(1:numberElements-1);
   phi(:,n+1)=x_sol(numberElements:end); 
end
w1=zeros(numberElements+1,M+1);
w1(activeDof_w,:)=w;
phi1=zeros(numberElements+1,M+1);
phi1(activeDof_phi,:)=phi;
end