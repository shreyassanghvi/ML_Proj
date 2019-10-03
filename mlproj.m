%mlproj(n,m) = mlproj(dimension , no of data points)
%
function mlproj(n,m,h)
    if nargin <3
        h=1;
    end
    X=rand(m,n);%use randn
    Avg=(sum(sum(X)))/(size(X,1)*size(X,2));
    K=(Avg-X)/h
    f=K/(m*(h^n))
    min(f)
end 