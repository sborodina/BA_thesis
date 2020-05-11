function u=boundary(u,ja,jb,ib)
 if ib==1        % flat BC
     u(1:ja-1,:)=ones(ja-1,1)*u(ja,:);
     u(jb+1:jb+ja-1,:)=ones(ja-1,1)*u(jb,:); 
 end
end