function NCMHap_poly(varargin)
% Implementation of NCMHap triploid(Zamani 2019)
clc
clear
if nargin==0
    Mpath=uigetdir('F:\znu\thesis\thesis\code\DataSet\Althap-Dataset\Data_100\P3-A2\sim-cov10-err0.1\','Select Dataset'); %one can change the file path and file name(for one length value).
    directory=[Mpath,'\sim*.txt'];
    D1=dir(directory);%read whole file and folder
    if isempty(D1)
        fprintf(['ERROR! Fail to find ',Mpath,'\n']);
        return;
    end
    %read real hap
    Hpath=uigetdir('F:\znu\thesis\thesis\code\DataSet\Althap-Dataset\Data_100\P3-A2\sim-cov10-err0.1\','Select Real Haplotype'); %one can change the file path and file name(for one length value).
    Hdirectory=[Hpath,'\sim*-phased.txt'];
    HD1=dir(Hdirectory);%read whole file and folder
    if isempty(HD1)
        fprintf(['ERROR! Fail to find ',Hpath,'\n']);
        return;
    end

    ploidy = 3;
    maxit = 100;    % Max. iteration
    expo = 2;		% Exponent for U		
    display = 1;	% Display info or not
    tStart=tic;
    cluster_n = ploidy;
    min_impro = 1e-5;		% Min. improvement
end
   
Num_paras=length(D1);
rrs=zeros(1,Num_paras);
finalmec=zeros(1,Num_paras);
MEC_n=1;
for d=1:Num_paras
    D1(d).name%for each parameter set
    %ss=[Mpath,D1(d).name,'\*.m.4.err']
    D=dir([Mpath,'\',D1(d).name]);
    if isempty(D)
        fprintf('ERROR! Fail to find %s *.txt!\n',D1(d).name);
        return;
    end
    HD=dir([Hpath,'\',HD1(d).name]);%read realhap
    if isempty(HD)
        fprintf('ERROR! Fail to find %s *.txt!\n',D1(d).name);
        return;
    end
    fprintf('Succeed to find %s!\n',D1(d).name);
    filename=[Mpath,'\',D1(d).name];
    fidin=fopen(filename);
    if fidin==-1
        fprintf('ERROR£¡fail to open the M file %d!\n',k);
        return;
    end
    Hfilename=[Hpath,'\',HD1(d).name];
    Hfidin=fopen(Hfilename);
    if Hfidin==-1
        fprintf('ERROR£¡fail to open the M file %s!\n',HD1(d).name);
        return;
    end
%     Hreal = textscan(Hfidin,'%s','delimiter','\n');
%     H_real=char(Hreal{1});
%     fclose(Hfidin);
    
    %% Reading the data
    start = 3;
    file = textscan(fidin,'%s','delimiter','\n');
    content = char(file{1,1});
    [read_num, ~] = size(content);
    read_num = read_num-2;
    line = textscan(content(2,:),'%s','delimiter',' ');
    hap_len = str2double(line{1}{1});
    frags=repmat('-',read_num,hap_len);
    for i = start:read_num+2
        line = textscan(content(i,:),'%s','delimiter',' ');
        line = line{1};
        for j = 1:str2double(line{1})
            frags(i - start + 1,str2double(line{j*2+1})+1:str2double(line{j*2+1})+length(line{j*2+2})) = ...
                line{j*2+2};
        end
    end
    
%     fileID = fopen(Hfidin,'r');
    sizeA = [3 Inf];
    Hreal = fscanf(Hfidin,'%s',sizeA);
    
    
    alleles_num = length(unique(frags))-1;
    [read_num, read_len] = size(frags);
    R = char();
    for i=1:read_num
        for j=1:read_len
            if frags(i,j)=='1'
                R(i,j)='a';
            elseif frags(i,j)=='0'
                R(i,j)='t';
            else
                R(i,j)='-';
            end
        end
    end
    
    [hap_num, hap_len] = size(Hreal);
    H_real = char();
    for i=1:hap_num
        for j=1:hap_len
            if Hreal(i,j)=='1'
                H_real(i,j)='a';
            elseif Hreal(i,j)=='0'
                H_real(i,j)='t';
            else
                H_real(i,j)='-';
            end
        end
    end
    
    m=size(R,1);
    BestChr=zeros(1,m);
%% GAHAP
        PS=500;
        j=1;
        MNG=1000;
        CHR=Create_rand(PS,m);
        previous=2000;
        thr=0;
        K=30;
        Best=1;
        Pc=0.75; % Crossover rate
        Pm=0.01; % Mutation rate
        while (j<=MNG && thr<=K && Best>0)
            if mod(j,50)==0
                fprintf('\n iteration:%d\t best=%d',j,Best);
            end
            F=HapFitness(CHR,frags);
            Best=min(F);
            index=find(F==Best);
            BestChr=CHR(index(1),:);            
            MatingPool=Rank(F);
            NewGen=SinglePoint(MatingPool,CHR,Pc);
            CHR=mutation(NewGen,Pm);
            j=j+1;
            if Best==previous
                thr=thr+1;
            else
                thr=0;
                previous=Best;
            end
        end
    %% create initial Haplotype    
    H_rec=MakeHapbyMajority(BestChr,R,cluster_n);
    H_rec=ImproveHap(H_rec,BestChr,R,cluster_n);
    %% NCMHap
    [T,I,F] = initNCMpoly(H_rec,R,cluster_n);% Initial fuzzy partition
    data=mapfag(R);
    delta=25;
    % Main loop
    for j = 1:maxit,
        [T,I,F, center, obj_fcn(j)] = stepNCMhap(data, T,I,F,cluster_n,delta,expo);
        if display,
            fprintf('Iteration count = %d, obj. fcn = %f\n', j, obj_fcn(j));
        end
        % check termination condition
        if j > 1,
            if abs(obj_fcn(j) - obj_fcn(j-1)) < min_impro, break; end,
        end
    end

    iter_n = j;	% Actual number of iterations
    [H_rec]=maphappoly(center);%return number  to allel
    rrs(MEC_n)=RR_ploidy(H_real,H_rec);% calculate Recostruction Rate
    finalmec(MEC_n)=PloidMEC(H_rec,R);% calculate MEC Rate
    tElapsed(MEC_n)=toc(tStart);% calculate time
    fprintf('The instance %s is finished!\n',D1(d).name);
    fid=fopen(['NCM_',D1(d).name,'_RR.txt'],'w+');
    fid2=fopen(['NCM_',D1(d).name,'_MEC.txt'],'w+');

    fprintf(fid,'%d\t%1.4f\t%d\t\n',j,rrs(MEC_n),tElapsed(MEC_n));
    fprintf(fid2,'%d\t%1.4f\t%d\t\n',j,finalmec(MEC_n),tElapsed(MEC_n));
    fclose(fid);
    fclose(fid2);
    MEC_n=MEC_n+1;
end
fid=fopen('NCM_all_RR.txt','w+');
fid2=fopen('NCM_all_MEC.txt','w+');
average_rr=mean(rrs);
average_mec=mean(finalmec);
average_t=mean(tElapsed);
fprintf(fid,'%1.4f\t%1.4f\t%1.4f\n',rrs,average_rr,average_t);
fprintf(fid2,'%1.4f\t%1.4f\t%1.4f\n',finalmec,average_mec,average_t);
fclose(fid);
fclose(fid2);

end
