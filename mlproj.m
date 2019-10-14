%mlproj(n,m) = mlproj(mean 1, varience 1, mean 2, varience 2, samples 1, samples 2, dimension, smoothing parameter)
%
function mlproj(u1,v1,u2,v2,p1,p2,n,h)
    if nargin <7
        h=1;
    end
    r=0.05;
    mr=(p1+p2)*r;
    I=[0.5+mr];
    I=round(I);
    r1=u1+v1.*randn(p1,n);
    r2=u2+v2.*randn(p2,n);
    X=[r1 ; r2 ];%use randn
    %X=[1 2;10 12; 7 11];
    Fx=ones(size(X));
    for i=1:(size(X,1))
        tmp=[X(1:i-1,:);X(i+1:end,:)];
        
        if (i==size(X,1))
            tmp=[X(1:i-1,:)];
        end
        t=X(i,:)-tmp;
        K=exp((t.*t)/-2)/sqrt(2*pi);
        Fx(i,:)=prod(sum(K),2);
    end
    Fx=Fx(:,1);
    Func=Fx/((p1+p2)*(h^(n)));
    Func_sort=sort(Func);
    if(mr<0.5)
        q=Func_sort(1)
    else
        q=(0.5+I-mr)*Func_sort(I)+(0.5-I+mr)*Func_sort(I+1)
    end
    Z=Func<=q
end 