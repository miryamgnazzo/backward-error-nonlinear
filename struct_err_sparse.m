function [D0,D1,D2]=struct_err_sparse(A0,A1,A2,V,L,k)

%initial step
[E0,E1,E2]=unstruct_error(A0,A1,A2,V,L);
D0=(A0 ~= 0) .* (E0+A0);
D2=(A2 ~= 0) .* (E2+A2);
D1=(A1 ~= 0) .* (E1+A1);

for j=1:50
    [E0,E1,E2]=unstruct_error(D0,D1,D2,V,L);
    D0=(A0 ~= 0) .* (E0+A0);
    D2=(A2 ~= 0) .* (E2+A2);
    D1=(A1 ~= 0) .* (E1+A1);
    
    fprintf('iteration %d %e \n',j,norm([E0 E1 E2],'fro'))
end