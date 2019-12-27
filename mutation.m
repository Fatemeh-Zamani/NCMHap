function NewGen=mutation(G,Pm)
[R,C]=size(G);

for r=1:R
    for c=1:C
        p=rand(1);
        if p<Pm
            if G(r,c)==0
                G(r,c)=1;
            else
                G(r,c)=0;
            end
        end
    end
end
NewGen=G;

end

