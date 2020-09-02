function U = initdist(Rh1,Rh2,data)
[R]=size(data,1);%col size- number of SNP
for r=1:R
    U(1,r)=FuzzyDis(Rh1,data(r,:));
    U(2,r)=FuzzyDis(Rh2,data(r,:));
end

