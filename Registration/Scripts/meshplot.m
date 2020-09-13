

function meshplot(x,y,n,m)

if ~exist('n','var') || isempty(n)
    [m, n] = size(x); 
end;

xy=[x(:) y(:)];
A = spdiags(ones(m*n,2),[1 m],m*n,m*n);
a=m:m:size(xy,1)-1;
A(sub2ind([m*n m*n],a,a+1))=0;
gplot(A,xy);axis ij equal off;
