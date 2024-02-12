function [D0,D1,D2]=unstruct_error(A0,A1,A2,V,L)
%calcolo delle perturbazioni 

p=size(V,2);
n=size(A0,1);

Ld=diag(L);

R=A2*V*L^2+A1*V*L +A0*V;
%r=R(:);

F=zeros(p,3);

for j=1:3
    F(:,j)=Ld.^(j-1);
end

K=Khatri_Rao(F,V');

delta=-R*pinv(K');
%keyboard

D0=delta(:,1:n);
D1=delta(:,n+1:2*n);
D2=delta(:,2*n+1:3*n);

