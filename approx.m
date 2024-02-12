function [V,L,index]=approx(A0,A1,A2)
%find approximation of the eigenpairs of the matrix polynomials

[VV,LL]=polyeig(A0,A1,A2);


[nn,I]=sort(abs(imag(((LL)))));

LL=LL(I);
VV=VV(:,I);

nr= sum(imag(LL)==0);

index=[1:nr , nr+2:2:size(VV,2)];

V=zeros(size(VV));
L=zeros(length(LL),length(LL));

V(:,1:nr)=VV(:,1:nr);
L(1:nr,1:nr)=diag(LL(1:nr));

j=nr+1;

while (j<size(VV,2))
   Vj=[real(VV(:,j))  -imag(VV(:,j))];
   Lj=[real(LL(j)) -imag(LL(j)); imag(LL(j)) real(LL(j))]; 
   
   V(:,j:j+1)=Vj;
   L(j:j+1,j:j+1)=Lj;
   
   j=j+2;
    
end

