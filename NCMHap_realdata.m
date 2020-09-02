function [rrs,MECscore] = NCMHap_realdata(varargin)
% Implementation of NCMHap(Zamani 2019)
% tStart=tic;     %calculate CPU time
if nargin==4
    M0=varargin{1};
    H_real=varargin{2};
    i=varargin{3};
    [m0,n0]=size(M0);
    Original=varargin{4};
end

Pivot1=1;
expo = 2;		% Exponent for U
max_iter = 50;		% Max. iteration
min_impro = 1e-5;		% Min. improvement
display = 1;		% Display info or not
rrs=[];
MECscore=[];
cluster_n=2;

[M,m,n,H_assem,Most]=Heter4to2(M0,m0,n0);
if n0<=2
    rrs=1;
    MECscore=0;%calculate Recostruction Rate
    return
end

% [M,m,n,H_assem,Most]=Heter4to2(M0,m0,n0);
if n==0
    rrs=1;
    MECscore=0;%calculate Recostruction Rate
    return
end
BestChr=zeros(1,m);
[C1,C2,~]=FCGraph(M,Pivot1);
BestChr(C2)=1;

H_rec=MakeHapbyMajority(BestChr,M);
H_rec=ImproveHap(H_rec,BestChr,M);

H_rec2=H_rec;
for ind=1:length(H_rec)
    if H_rec(ind)=='t'
        H_rec2(ind)='a';
    else
       H_rec2(ind)='t';
    end
end
H_rec=[H_rec;H_rec2];

obj_fcn = zeros(max_iter, 1);
data=mapfag(M);


[T,I,F] = initNCM(H_rec,M,cluster_n);% Initial fuzzy partition
delta=25;
% Main loop
for j = 1:max_iter,
    [T,I,F, center, obj_fcn(j)] = stepNCM(data, T,I,F,cluster_n,delta,expo);
    if display,
        fprintf('Iteration count = %d, obj. fcn = %f\n', j, obj_fcn(j));
    end
    % check termination condition
    if j > 1,
        if abs(obj_fcn(j) - obj_fcn(j-1)) < min_impro, break; end,
    end
end

% tElapsed=toc(tStart);% calculate CPU time
[H_rec,H_rec2]=maphap(center);%return number  to allel
H_assem1=CompleteHap(H_rec,H_assem);
H_assem4=H_assem2to4(H_assem1,Most);
rrs=RR4(H_real,H_assem4);% calculate Recostruction Rate
MECscore=MEC(H_assem4(1,:),M0);%calculate Recostruction Rate
%end for each parameter set

fprintf('one paramenter set of one length value ');
end
