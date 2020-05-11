% Nessyahu-Tadmor scheme for BBM equation
% u_t+uu_x-u_xxt=0
% K_t+(.5*u^2)_x=0 
% K=u-u_xx   

clear
global dx 

N=6000;          % Number of grid points
xa=-210; xb=70;  % Computational domain
tmax=300;        % Final computational time 
ic=1;            % Initial condition  
ib=1;            % Boundary condition 
islope=1;        % Slope limiter 
CFL=.4;          % CFL condition

% ========================================================
ja=2+islope; jb=N+ja-1; Ntot=jb+ja-1;
J=ja:jb; JJ=ja-1:jb+1; Jm=ja-1:jb;
dx=(xb-xa)/N; x=xa+dx/2:dx:xb-dx/2; % Grid

t=0; u0=initial(x,ic); % Initial condition
vel0=u0(:,1);
K0=u0(:,2);

figure
hold on
title(['Time   ',num2str(t)]);
plot(x,vel0);
xlabel('x'); ylabel('u');
hold off

u(J,1)=K0; 
vel(J,1)=vel0; 

[nn, kk]=size(K0);
u1=zeros(Ntot,kk); f1=u1; um=u1; 
fm=u1; w=u1; velm=u1;

tp=0; dtp=.3;

% ========================================================
 while (t<tmax)
   dt=CFL/maxeig(vel(J))*dx;
   if (t+2*dt>tmax)
      dt=(tmax-t)/2;
   end
   lambda=dt/dx;
  
   u=boundary(u,ja,jb,ib);
   vel=boundary(vel,ja,jb,ib);
   u1(JJ,:)=slope(u,ja-1,jb+1,islope);
   f=flux(vel); 
   f1(JJ,:)=slope(f,ja-1,jb+1,islope);
   um(JJ,:)=u(JJ,:)-lambda*f1(JJ,:)/2;  
   vel=boundary(vel,ja,jb,ib);
   velm(JJ,:)=sweepp(um(JJ,:),vel(ja-2,:),vel(ja-1,:));  
   fm(JJ,:)=flux(velm(JJ,:));
   w(Jm,:)=(u(Jm+1,:)+u(Jm,:))/2+(u1(Jm,:)-u1(Jm+1,:))/8 ...
                 -lambda*(fm(Jm+1,:)-fm(Jm,:));   
   velm=boundary(velm,ja,jb,ib);            
   vel(Jm,:)=sweepp(w(Jm,:),velm(ja-2,:),velm(ja-1,:)); 
                                
   w=boundary(w,ja,jb,ib);
   vel=boundary(vel,ja,jb,ib);
   u1(Jm,:)=slope(w,ja-1,jb,islope);
   f=flux(vel);
   f1(Jm,:)=slope(f,ja-1,jb,islope);
   um(Jm,:)=w(Jm,:)-lambda*f1(Jm,:)/2;
   vel=boundary(vel,ja,jb,ib);
   velm(Jm,:)=sweepp(um(Jm,:),vel(ja-2,:),vel(ja-1,:));
   fm(Jm,:)=flux(velm(Jm,:));
   u(J,:)=(w(J-1,:)+w(J,:))/2-(u1(J,:)-u1(J-1,:))/8 ...
           -lambda*(fm(J,:)-fm(J-1,:));
   vel(J,:)=sweepp(u(J,:),velm(ja-1,:),velm(ja,:));

   if t>tp
       %plot(x,vel(J,1),x,u(J,1));
       plot(x,vel(J,1));
       axis([xa xb -2 .2]);
       title(['Time   ',num2str(t)])
       pause(0.01)
       tp=tp+dtp;
   end
   % ============================== 
   t=t+2*dt;   
 end
% ========================================================

figure
hold on
title(['Time   ',num2str(t)]);
plot(x,vel(J,1));
axis([xa xb -2 .2]);
xlabel('x'); ylabel('u');
legend('u_L=0, u_R=-1');
hold off





