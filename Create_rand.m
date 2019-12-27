function  CHR=Create_rand(PS,L)
CHR=zeros(PS,L);
% for i=1:PS
%     for j=1:L
%         CHR(i,j)=randi([0 2],5,2);
%     end
% end
CHR=randi([0 2],PS,L);
end

