function s = slope(u,ja,jb,is)
if is==1 % Simple MinMod limiter
   J = ja:jb;
   s = mm(u(J+1,:)-u(J,:),u(J,:)-u(J-1,:));
end
if is==2 % UNO limiter
   J = ja:jb;
   JJ = ja-1:jb+1;
   Jm = ja-1:jb;
   D(JJ,:) = u(JJ+1,:)-2*u(JJ,:)+u(JJ-1,:);
   d(Jm,:) = u(Jm+1,:)-u(Jm,:);
   s = mm(d(J-1,:)+0.5*mm(D(J-1,:),D(J,:)),d(J,:)-0.5*mm(D(J,:),D(J+1,:)));
end
