function score = PloidMEC(Haps,M)
%#codegen
score=0;
[R,~]=size(M);
n=0;
[Nhap,~]=size(Haps);
d=zeros(1,Nhap);
for r=1:R
    for i=1:Nhap
        d(i)=dist_of_seq(Haps(i,:),M(r,:));
    end
    score=score+min(d);
end
end

