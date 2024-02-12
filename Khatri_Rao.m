function K=Khatri_Rao(F,V)

K=zeros(size(F,1),size(F,2)*size(V,2));

if (size(F,1)~=size(V,1))
    error('inconsistent dimension') 
end

for i=1:size(F,1)
    K(i,:)=kron(F(i,:),V(i,:)); 
end
