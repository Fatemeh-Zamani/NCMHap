function [Rh]=maphappoly(center)
R=size(center,2);
for r=1:R
    if center(1,r)>0
        Rh1(r)='a';
    else
        Rh1(r)='t';
    end
    if center(2,r)>0
        Rh2(r)='a';
    else
        Rh2(r)='t';
    end
    if center(3,r)>0
        Rh3(r)='a';
    else
        Rh3(r)='t';
    end
end
Rh=[Rh1;Rh2;Rh3];
end

