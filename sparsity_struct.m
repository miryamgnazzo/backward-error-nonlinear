function P=sparsity_struct(A)
%it reconstructs the matrix P associated with the matrix A
%consider as structure the sparsity pattern induced by A

[m,n]=size(A);
[row,column,~]=find(A);

s=length(row);

P=zeros(m*n,s);

for i = 1 : s
   
    index=(column(i)-1)*m + row(i);

    P(index,i) = 1; 
end