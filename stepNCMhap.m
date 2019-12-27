function [T_new, I_new,F_new, center, obj_fcn] = stepNCMhap(data, T,I,F,cluster_n,delta, expo)
% mf = U.^expo;       % MF matrix after exponential modification
mf = T.^expo  ;   
Iexpo = I.^expo;
Fexpo =F.^expo;

data_num = size(data,1);
data_dim = size(data,2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%new center (obtained haplotype)
 center=makecenter(data,T,cluster_n);
% center = mf*data./(sum(mf,2)*ones(1,size(data,2))); %new center
dist = distfcm(center, data);       % fill the distance matrix
obj_fcn1 = sum(sum((dist.^2).*mf));  % objective function
% tmp = dist.^(-2/(expo-1));      % calculate new U, suppose expo != 1
% U_new = tmp./(ones(2, 1)*sum(tmp));
if cluster_n==2
    mvec=mean(center,1);
    dist1 = distfcm(mvec, data);
    obj_fcn2 = (sum((dist1.^2).*(Iexpo)));  
    obj_fcn3 = sum(delta^2.*Fexpo);
    obj_fcn=obj_fcn1+obj_fcn2+obj_fcn3;
      K_tmp = (.7*sum(dist.^(-2/(expo-1)))+.2*(dist1.^(-2/(expo-1)))+.1*delta^(-2/(expo-1)));
    T_new = .7*(dist).^(-2/(expo-1))./(ones(cluster_n, 1)*K_tmp);
    I_new = .2*(dist1).^(-2/(expo-1))./K_tmp;
    F_new = .1*(delta).^(-2/(expo-1))./K_tmp;


else
    tg=sort(T,'descend');
    tg=tg(1:2,:);
    for i=1:data_num
        for j=1:2
            a=find(T(:,i)==tg(j,i));
            if tg(1,i)==tg(2,i)
                if length(a)==2
                    c_ind(j,i)=a(1);
                    c_ind(j+1,i)=a(2);
                    break
                elseif length(a)==3
                    c_ind(j,i)=a(1);
                    c_ind(j+1,i)=a(2);
                    c_ind(j+2,i)=a(3);
                    break
                end
            end
            if length(a)==2
                c_ind(j,i)=a(1);
                c_ind(j+1,i)=a(2);
%                 break
            elseif length(a)==3
                c_ind(j,i)=a(1);
                c_ind(j+1,i)=a(2);
                c_ind(j+2,i)=a(3);
                break
            else
                c_ind(j,i)=a;
            end
        end
        c_ind2=zeros();
        for kkk=1:size(c_ind(:,i))
            if c_ind(kkk,i)>0
                c_ind2(kkk,1)=c_ind(kkk,i);
            end
        end
%         i
%         c_ind2
        cen_mn(i,:)=mean([center(c_ind2(:,1),:)]);    
    end
    
    dist1=sqrt(sum(((data-cen_mn).^2)'));
    obj_fcn2 = (sum((dist1.^2).*(Iexpo)));
    obj_fcn3 = sum(delta^2.*Fexpo);
    obj_fcn=obj_fcn1+obj_fcn2+obj_fcn3;
     K_tmp = (.7*sum(dist.^(-2/(expo-1)))+.2*(dist1.^(-2/(expo-1)))+.1*delta^(-2/(expo-1)));
    T_new = .7*(dist).^(-2/(expo-1))./(ones(cluster_n, 1)*K_tmp);
    I_new = .2*(dist1).^(-2/(expo-1))./K_tmp;
    F_new = .1*(delta).^(-2/(expo-1))./K_tmp;

end