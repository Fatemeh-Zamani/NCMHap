function center=makecenter(data,T,cluster_n)
if cluster_n==2
    C=size(data,2);
    maxT = max(T);

    index1 = find(T(1,:) == maxT);
    index2 = find(T(2,:) == maxT);
    center=zeros(2,C);
    for c=1:C
        center(1,c)=T(1,index1)*data(index1,c);
        center(2,c)=T(2,index2)*data(index2,c);

        center(1,c)=center(1,c)/(length(index1)+1);
        center(2,c)=center(2,c)/(length(index2)+1);
        if center(1,c)>0
            center(1,c)=1;
        else
            center(1,c)=-1;
        end
        if center(2,c)>0
            center(2,c)=1;
        else
            center(2,c)=-1;
        end
    end
elseif cluster_n==3
    C=size(data,2);
    maxT = max(T);

    index1 = find(T(1,:) == maxT);
    index2 = find(T(2,:) == maxT);
    index3 = find(T(3,:) == maxT);
    center=zeros(3,C);
    for c=1:C
        center(1,c)=T(1,index1)*data(index1,c);
        center(2,c)=T(2,index2)*data(index2,c);
        center(3,c)=T(3,index3)*data(index3,c);

        center(1,c)=center(1,c)/(length(index1)+1);
        center(2,c)=center(2,c)/(length(index2)+1);
        center(3,c)=center(3,c)/(length(index3)+1);
        if center(1,c)>0
            center(1,c)=1;
        else
            center(1,c)=-1;
        end
        if center(2,c)>0
            center(2,c)=1;
        else
            center(2,c)=-1;
        end
        if center(3,c)>0
            center(3,c)=1;
        else
            center(3,c)=-1;
        end
    end
end


end

