function data=mapfag(frag)
[R,C]=size(frag);
for r=1:R
    for c=1:C
        if frag(r,c)=='a'
            data(r,c)=1;
        elseif frag(r,c)=='t'
            data(r,c)=-1;
        else
            data(r,c)=0;
        end
    end
end
end

