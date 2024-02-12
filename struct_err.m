function [D0,D1,D2]=struct_err(A0,A1,A2,V,L,k)

%initial step
[E0,E1,E2]=unstruct_error(A0,A1,A2,V,L);
D0=E0+A0;
D2=E2+A2;
D1=E1+A1;
[UU,SS,VV]=svd(D1);
D1=UU(:,1:k)*SS(1:k,1:k)*VV(:,1:k)';

for j=1:50
    [E0,E1,E2]=unstruct_error(D0,D1,D2,V,L);
    D0=E0+D0;
    D2=E2+D2;
    D1=E1+D1;
    [UU,SS,VV]=svd(D1);
    D1=UU(:,1:k)*SS(1:k,1:k)*VV(:,1:k)';
    
    fprintf('iteration %d %e \n',j,norm([E0 E1 E2],'fro'))
end