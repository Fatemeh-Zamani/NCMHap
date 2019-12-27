function  NewCHR=SinglePoint(Parents,CHR,Pc)
[R,C]=size(CHR);
NewCHR=zeros(R,C);
N=length(Parents)*2;

for n=1:2:N
    p1=0;
    while p1==0
        p1=round(rand(1)*N/2);
    end
    p2=0;
    while p2==0 || p1==p2
        p2=round(rand(1)*N/2);
    end
    p=rand(1);
    parent1=CHR(Parents(p1),:);
    parent2=CHR(Parents(p2),:);
    if p<=Pc
        point=0;
        while point==0
        point=round(rand(1)*C);
        end
        c1(1:point)=parent1(1:point);
        c1(point+1:C)=parent2(point+1:C);
        c2(1:point)=parent2(1:point);
        c2(point+1:C)=parent1(point+1:C);
    else
        c1=parent1;
        c2=parent2;
    end
    NewCHR(n,:)=c1;
    NewCHR(n+1,:)=c2;
    

end

