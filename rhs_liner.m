function F=rhs_liner(q,N,L,activeDof_phi)
x_nodes=0:L/N:L;
b=zeros(N+1,1);
h=L/N;
% for i=1:N %N
%     xa=x_nodes(i);
%     xb=x_nodes(i+1);
%     xm=(xa+xb)/2;
%     f1=@(x)q(xm+h*x/2).*(1-x)/2;
%     f2=@(x)q(xm+h*x/2).*(1+x)/2;
%     be(1,1)=h*quadgk(f1,-1,1)/2;
%     be(2,1)=h*quadgk(f2,-1,1)/2;
%     b(i:i+1,1)=b(i:i+1,1)+be;
% end
cnt=N/2+1;
b(cnt)=1;
% F=b;
F=b(activeDof_phi);
end