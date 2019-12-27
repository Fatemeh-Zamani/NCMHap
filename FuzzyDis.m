function out=FuzzyDis(str1,str2)
% called by FuzzyPartion
L=length(str1);
d=0;m=0;
for l=1:L
    if str1(l)~=str2(l) && str1(l)~='-' && str2(l)~='-'
        d=d+1;
        m=m+1;
    elseif str1(l)=='-' && str2(l)~='-' || str2(l)=='-' && str1(l)~='-'
        d=d+0.5;
        m=m+1;
    elseif str1(l)==str2(l) && str1(l)~='-' && str2(l)~='-'
        m=m+1;
    end
end
out=d/m;

