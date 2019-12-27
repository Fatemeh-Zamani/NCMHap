function Hap=MakeHapbyMajority(BestChr,Frags,cluster_n)
if cluster_n==2
    [R,C]=size(Frags);
    [N,L]=size(BestChr);
    for n=1:N
        n1=find(BestChr(n,:)==0);
        n2=find(BestChr(n,:)==1);
        c1='';
        c2='';
        c1=Frags(n1,:);
        c2=Frags(n2,:);
        for cc=1:C
            a1=find(c1(:,cc)=='a');
            t1=find(c1(:,cc)=='t');
            a2=find(c2(:,cc)=='a');
            t2=find(c2(:,cc)=='t');
            if length(a1)+length(t2)>length(t1)+length(a2)
                Hap(cc)='a';
            else
                Hap(cc)='t';
            end
        end
    end
elseif cluster_n==3
    [R,C]=size(Frags);
    [N,L]=size(BestChr);
    for n=1:N
        n1=find(BestChr(n,:)==0);
        n2=find(BestChr(n,:)==1);
        n3=find(BestChr(n,:)==2);
        c1='';
        c2='';
        c3='';
        c1=Frags(n1,:);
        c2=Frags(n2,:);
        c3=Frags(n3,:);
        for cc=1:C
            a1=find(c1(:,cc)=='a');
            t1=find(c1(:,cc)=='t');
            a2=find(c2(:,cc)=='a');
            t2=find(c2(:,cc)=='t');
            a3=find(c3(:,cc)=='a');
            t3=find(c3(:,cc)=='t');
            if length(a1)>length(t1)
                Hap1(cc)='a';
            else
                Hap1(cc)='t';
            end
            if length(a2)>length(t2)
                Hap2(cc)='a';
            else
                Hap2(cc)='t';
            end
            if length(a3)>length(t3)
                Hap3(cc)='a';
            else
                Hap3(cc)='t';
            end
        end
    end 
    Hap=[Hap1;Hap2;Hap3];
end

