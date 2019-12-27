function NCMHap(varargin)
% Implementation of NCMHap(Zamani 2019)
clc
clear
Pivot1=0.65;
expo = 2;		% Exponent for U
max_iter = 200;		% Max. iteration
min_impro = 1e-5;		% Min. improvement
display = 1;		% Display info or not
if nargin==0
    Mpath=uigetdir('E:\karshenasi\term8\haplotype\Data\SurveyDataset\SurveyDataset','Select Dataset'); %one can change the file path and file name(for one length value).
    directory=[Mpath,'\exp-100-*-h0.4'];
    D1=dir(directory);%read whole file and folder
end
if isempty(D1)
    fprintf(['ERROR! Fail to find ',Mpath,'\n']);
    return;
end
Num_paras=length(D1);
allavgRR=[];
allavgMEC=[];
for d=1:Num_paras
        D1(d).name%for each parameter set
        %ss=[Mpath,D1(d).name,'\*.m.4.err']
        D=dir([Mpath,'\',D1(d).name,'\*.m.4.err']);
        if isempty(D)
            fprintf('ERROR! Fail to find %s *.m.4.err!\n',D1(d).name);
            return;
        end
        fprintf('Succeed to find %s \\*.m.4.err!\n',D1(d).name);
        %Num_exps=length(D);
        Num_exps=100;% number of instance to read
        rrs=zeros(1,Num_exps);
        cluster_n=2;
        for i=1:Num_exps  %for each instance
            tStart=tic;
            obj_fcn = zeros(max_iter, 1);
            [M0,m0,n0,H_real]=inputdata(Mpath,D1(d).name,D(i).name,i); % set frags on M0 and target haps in H_real
            [M,m,n,H_assem,Most]=Heter4to2(M0,m0,n0);
            BestChr=zeros(1,m);
            [C1,C2,~]=FCGraph(M,Pivot1);
            BestChr(C2)=1;
              
            H_rec=MakeHapbyMajority(BestChr,M,cluster_n);
            H_rec=ImproveHap(H_rec,BestChr,M,cluster_n);
            
            if cluster_n==2
                H_rec2=H_rec;
                for ind=1:length(H_rec)
                    if H_rec(ind)=='t'
                        H_rec2(ind)='a';
                    else
                       H_rec2(ind)='t';
                    end
                end
                H_rec=[H_rec;H_rec2];
            end
            [T,I,F] = initNCM(H_rec,M,cluster_n);% Initial fuzzy partition
            data=mapfag(M);
            delta=25;
            % Main loop
            for j = 1:max_iter,
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
            [H_rec,H_rec2]=maphap(center);%return number  to allel
            H_assem=CompleteHap(H_rec,H_assem);%create second haplotype
            H_assem4=H_assem2to4(H_assem,Most);
            rrs(i)=RR4(H_real,H_assem4);% calculate Recostruction Rate
            tElapsed(i)=toc(tStart);% calculate time
            fprintf('The instance %d is finished!\n',i);
        end %end for each instance
        
        fid=fopen(['NCM_',D1(d).name,'_RR.txt'],'w+');
        
        for j=1:Num_exps
            fprintf(fid,'%d\t%1.4f\t%d\t\n',j,rrs(j),tElapsed(j));
        end
        
        average_rr=mean(rrs);
        fprintf(fid,'%1.4f\n\n',average_rr);
        average_t=mean(tElapsed);
        fprintf(fid,'%1.4f\t%1.4f\n',average_rr,average_t);
        fclose(fid);
end %end for each parameter set

fprintf('one paramenter set of one length value ends£¡');









