function [Matrix2,Matrix1]=Transform2Geraci(File_SNPs)
fid_SNP=fopen(File_SNPs,'r');%
Matrix2='';
Matrix1='';
counter=1;
flag=1;
index=0;
while ~feof(fid_SNP)
    tline=fgetl(fid_SNP);
    str=strsplit(tline);
    N=str2num(char(str(1)));
    if flag==1
        Start=str2num(char(str(3)));flag=0;
    end
    if N>2
        'new challenge'
    end
    if N==1
        pos1=str2num(char(str(3)))-Start+1;
        frag=char(str(4));
    elseif N==2 
        tmp='';
        frag1=char(str(4));
        pos1=str2num(char(str(3)))-Start+1;
        frag2=char(str(6));
        pos2=str2num(char(str(5)))-Start+1;
        tmp(1:pos2-(pos1+length(frag1)))='-';
        frag=[frag1,tmp,frag2];
    elseif N==3
        tmp='';
        tmp2='';
        frag1=char(str(4));
        pos1=str2num(char(str(3)))-Start+1;
        frag2=char(str(6));
        pos2=str2num(char(str(5)))-Start+1;
        tmp(1:pos2-(pos1+length(frag1)))='-';
        frag3=char(str(8));
        pos3=str2num(char(str(7)))-Start+1;
        tmp2(1:pos3-(pos2+length(frag2)))='-';
        frag=[frag1,tmp,frag2,tmp2,frag3];
    end
        M(counter)={frag};
        index(counter)=pos1;
        counter=counter+1;
end

N_frags=length(M);
C=index(N_frags)+length(char(M(N_frags)))-1;
Matrix2(1:N_frags,1:C)='-';
for i=1:N_frags
    df=index(i):index(i)+length(char(M(i)))-1;
    Matrix2(i,index(i):index(i)+length(char(M(i)))-1)=map2at(M(i));
end
[R,C]=size(Matrix2);
tmpStr(1:R)='-';
tmpMatrix=Matrix2';
i=1;
ind=[];
for c=1:C
    if tmpMatrix(c,:)==tmpStr
       ind(i)=c;
       i=i+1;
    end
end
%for j=1:length(ind)
    tmpMatrix(ind,:)=[];
%end
Matrix2=tmpMatrix';

fclose(fid_SNP);

