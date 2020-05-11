function m=mm(a,b)
m=min(a,b).*(a>0).*(b>0)+max(a,b).*(a<0).*(b<0);
