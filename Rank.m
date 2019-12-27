function Parents= Rank(fchr)
copy=fchr;
copy=sort(copy,'descend');
L=length(fchr);
for l=1:L
    for i=1:L
        if copy(l)==fchr(i)
            Rank(l)=i;
            fchr(i)=-1;
            break;
        end
    end
end
S=sum(Rank);
R=length(fchr);
Parents=zeros(1,round(R/2));
r=1;
while r<=round(R/2)
    p=rand(1);
    domain=0;
    i=1;
    while p>domain && i<=R
        domain=domain+i/S;
        i=i+1;
    end
    f=0;
    for j=1:r
        if Parents(j)==Rank(i-1)
            f=1;
        end
    end
    if f==0
        Parents(r)=Rank(i-1);
        r=r+1;
    end
    
    
end
