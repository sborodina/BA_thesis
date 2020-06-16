function vel=sweepp(u,v0,v1)
% u_xx-u=-K
% u_0=kappa1*u_1+mu1 
% u_N=kappa2*u_{N-1}+mu2
global dx
  N=length(u(:,1)); 
  kappa1=v0/v1; mu1=0; 
  kappa2=1; mu2=0;
  alpha=zeros(N,1); 
  betta=alpha;
  K=u(:,1); 
  a=ones(N,1); b=a;
  c=a+b+dx^2; f=K*dx^2;  
  alpha(1)=kappa1; betta(1)=mu1;
  for i=1:N-1
       d=c(i)-alpha(i)*a(i);
       alpha(i+1)=b(i)/d;
       betta(i+1)=(a(i)*betta(i)+f(i))/d;
  end
  vel(N)=(kappa2*betta(N)+mu2)/(1-kappa2*alpha(N));
  for i=1:N-1
       j=N-i;
       vel(j)=alpha(j+1)*vel(j+1)+betta(j+1);
  end
  vel=vel';
end
