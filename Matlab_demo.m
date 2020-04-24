%mlproj(n,m) = mlproj(mean 1, varience 1, mean 2, varience 2, samples 1, samples 2, dimension, smoothing parameter)
%
function ML(u1,v1,u2,v2,p1,p2,n,h)
    if nargin <7
        h=1;
    end
    m=p1+p2;
    hx=(8*sqrt(pi)*(1/(2*sqrt(pi))))/(3*m);
    r=0.05;
    mr=(p1+p2)*r;
    I=0.5+mr;
    I=floor(I);
    r1=u1+v1.*randn(p1,n);
    r2=u2+v2.*randn(p2,n);
    X=[r1 ; r2 ];%use randn
    sigma=sqrt(((1/(m-1))*sum(X.^2,1))-((1/(m*(m-1)))*(power(sum(X,1),2))));
    h=hx*sigma;
    Fx=ones(size(X));
    for i=1:(size(X,1))
        tmp=[X(1:i-1,:);X(i+1:end,:)];
        if (i==size(X,1))
            tmp=[X(1:i-1,:)];
        end
        t=((X(i,:)-tmp)./h);
        K=exp((t.*t)/-2)/sqrt(2*pi);
        Fx(i,:)=prod(sum(K),2);
    end
    Fx=Fx(:,1);
    Func=Fx/(m*(prod(h)));
    Func_sort=sort(Func);
    if(mr<0.5)
        q=Func_sort(1);
    else
        q=(0.5+I-mr)*Func_sort(I)+(0.5-I+mr)*Func_sort(I+1);
    end
    Z=Func<=q;
    si=find(Z,1);
    L = zeros(si,1);
    X1 = zeros(si,n);
    j=1;
    for i= 1:m
        if(Func(i)<=q)
            L(j) = Func(i); % function values of all Atypical elements stored in L
            X1(j) = X(i); % All point values of Atypical elements stored in X1
            j=j+1;
        end
    end
    si
    % L is the matrix with function values of the atypical elements
    %Extension of Population
    %1. Making the parameters for addition of more atypical elements
    a= min(X1 - ones(si,n))
    b= max(X1 + ones(si,n))
    c= max(L+ones(size(L,1),1),L-ones(size(L,1),1));
    c1 = max(c)
    L;
    typ = m - si; % Number of typical elements
    rem = typ- si; % Number of atypical elements to generate
    k=1;
    new_At = zeros(rem,n);
    while k <=rem
        u = a + (b-a).*rand(1,n);
        v= 0 + (c-0).*rand();
        t=((u-tmp)./h);
        K=exp((t.*t)/-2)/sqrt(2*pi);
        F_u=prod(sum(K),2);
        if F_u >= v
            for l=1:n
                new_At(k,l) = u(1,l);
            end
            k=k+1;
        end
    end
    All_At = [X1;new_At];
    All_At=[All_At, zeros(size(All_At,1),1)]; % labelling all atypical elements as 0
    r1=[r1 , ones(size(r1,1),1)]; % Labelling all typical elements as 1
    All = [All_At;r1]; %Complete Dataset with equal number of typical and atypical elements
    %Training
    X = All(1:size(All,1),1:n);
    Y= All(1:size(All,1),n+1);
    Mdl = fitctree(X,Y);
    view(Mdl);
    view(Mdl,'mode','graph','Prune',4)
end
