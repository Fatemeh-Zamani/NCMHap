function H_rec=ImproveHap(H_rec,BestChr,frags,cluster_n)
if cluster_n==2
    Hnew=H_rec;
    MecH_rec=MEC(H_rec,frags);
    Last=MecH_rec;
    MecHnew=0;
    while MecHnew<Last
        C1=find(BestChr==0);
        C2=find(BestChr==1);
        cluster1=frags(C1,:);
        cluster2=frags(C2,:);
        L=length(H_rec);
        % Find centers of clusters
        for cc=1:L
            a1=find(cluster1(:,cc)=='a');
            t1=find(cluster1(:,cc)=='t');
            a2=find(cluster2(:,cc)=='a');
            t2=find(cluster2(:,cc)=='t');
            if length(a1)==length(t1) || length(a2)==length(t2)
                disp('equal');
            end
            if length(a1)>length(t1)
                Center1(cc)='a';
            else
                Center1(cc)='t';
            end
            if length(a2)>length(t2)
                Center2(cc)='a';
            else
                Center2(cc)='t';
            end
        end
        % Find frags with max distances from their centers
        [R,~]=size(cluster1);
        dist1=zeros(1,R);
        for i=1:R
            dist1(i)=FuzzyDis(cluster1(i,:),Center1);
        end
        Maxdist1=max(dist1);
        f1=find(dist1==Maxdist1,1,'first');
        %+++++++++++++++++++++++++++++++++++++++++++
        [R,~]=size(cluster2);
        dist2=zeros(1,R);
        for i=1:R
            dist2(i)=FuzzyDis(cluster2(i,:),Center2);
        end
        Maxdist2=max(dist2);
        f2=find(dist2==Maxdist2,1,'first');

        BestChr(f1)=1;
        BestChr(f2)=0;

        Hnew=MakeHapbyMajority(BestChr,frags,cluster_n);
        MecHnew=MEC(Hnew,frags);
        if MecH_rec> MecHnew
            H_rec=Hnew;
            Last=MecH_rec;
            MecH_rec=MecHnew;
        else
            Last=MecHnew;
        end
    end
elseif cluster_n==3
    MecH_rec=PloidMEC(H_rec,frags);
    Last=MecH_rec;
    MecHnew=0;
    while MecHnew<Last
        C1=find(BestChr==0);
        C2=find(BestChr==1);
        C3=find(BestChr==2);
        cluster1=frags(C1,:);
        cluster2=frags(C2,:);
        cluster3=frags(C3,:);
        Hnew1=H_rec(1,:);
        L=length(Hnew1);
        % Find centers of clusters
        for cc=1:L
            a1=find(cluster1(:,cc)=='a');
            t1=find(cluster1(:,cc)=='t');
            a2=find(cluster2(:,cc)=='a');
            t2=find(cluster2(:,cc)=='t');
            a3=find(cluster3(:,cc)=='a');
            t3=find(cluster3(:,cc)=='t');
            if length(a1)>length(t1)
                Center1(cc)='a';
            else
                Center1(cc)='t';
            end
            if length(a2)>length(t2)
                Center2(cc)='a';
            else
                Center2(cc)='t';
            end
            if length(a3)>length(t3)
                Center3(cc)='a';
            else
                Center3(cc)='t';
            end
        end
        % Find frags with max distances from their centers
        [R,~]=size(cluster1);
        dist1=zeros(1,R);
        for i=1:R
            dist1(i)=FuzzyDis(cluster1(i,:),Center1);
        end
        Maxdist1=max(dist1);
        f1=find(dist1==Maxdist1,1,'first');
        %+++++++++++++++++++++++++++++++++++++++++++
        [R,~]=size(cluster2);
        dist2=zeros(1,R);
        for i=1:R
            dist2(i)=FuzzyDis(cluster2(i,:),Center2);
        end
        Maxdist2=max(dist2);
        f2=find(dist2==Maxdist2,1,'first');
        %+++++++++++++++++++++++++++++++++++++++++++
        [R,~]=size(cluster3);
        dist3=zeros(1,R);
        for i=1:R
            dist3(i)=FuzzyDis(cluster3(i,:),Center3);
        end
        Maxdist3=max(dist3);
        f3=find(dist3==Maxdist3,1,'first');

        distf1center2 = FuzzyDis(cluster1(f1,:),Center2);
        distf1center3 = FuzzyDis(cluster1(f1,:),Center3);
        if Maxdist1~=distf1center2 && Maxdist1~=distf1center3
            if distf1center2<distf1center3
                BestChr(f1)=1;
            else
                BestChr(f1)=2;
            end
        end
        
        distf2center1 = FuzzyDis(cluster2(f2,:),Center1);
        distf2center3 = FuzzyDis(cluster2(f2,:),Center3);
        if Maxdist2~=distf2center1 && Maxdist2~=distf2center3
            if distf2center1<distf2center3
                BestChr(f2)=0;
            else
                BestChr(f2)=2;
            end
        end
        distf3center1 = FuzzyDis(cluster3(f3,:),Center1);
        distf3center2 = FuzzyDis(cluster3(f3,:),Center2);
        if Maxdist3~=distf3center2 && Maxdist3~=distf3center1
            if distf3center1<distf3center2
                BestChr(f3)=0;
            else
                BestChr(f3)=1;
            end
        end
        % create new hap
        Hnew=MakeHapbyMajority(BestChr,frags,cluster_n);
        %calculate MEC for new hap
        MecHnew = PloidMEC(Hnew,frags);
        if MecH_rec> MecHnew
            H_rec=Hnew;
            Last=MecH_rec;
            MecH_rec=MecHnew;
        else
            Last=MecHnew;
        end
    end
end
end

