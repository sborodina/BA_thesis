function u=initial(x,k)
% initial velocity (vel) and variable (K)
 N=length(x); x0=0; 
 if k==1                  % Tipical Riemann problem
      uL=[-.1, -.1]; uR=[-1, -1]; nc=length(uL);
      u=reshape(kron(uL,(x<=x0)) + kron(uR,(x>x0)),N,nc);
 end
end
