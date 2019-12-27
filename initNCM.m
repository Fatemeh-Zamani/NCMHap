function [Tn,In,Fn,top] = initNCM(Rh,data,cluster_n)

if cluster_n==2
    Rh1=Rh(1,:);
    Rh2=Rh(2,:);
end
[R]=size(data,1);%col size- number of SNP
%for n cluster?
for r=1:R
    if cluster_n==2
        T(1,r)=FuzzyDis(Rh1,data(r,:));
        T(2,r)=FuzzyDis(Rh2,data(r,:));
    end
end     
% T = rand(cluster_n, data_n);
I=rand(1,R); 
F=rand(1,R);
% T=T.^(-2); 
T=[T;I;F];

% tmp = T.^(-2);  
col_sum = sum(T);
T = T./col_sum(ones(size(T,1), 1), :);
Tn=T(1:cluster_n,:);
In=T(cluster_n+1,:);
Fn=T(end,:);

