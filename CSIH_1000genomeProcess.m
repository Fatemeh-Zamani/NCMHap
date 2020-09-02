clc;clear;
current=pwd;
sizefrag=0;
prefix = 'F:\znu\thesis\thesis\code\DataSet\1000genome\1000genome\';
outprefix = 'F:\znu\thesis\thesis\code\DataSet\results\results\';

A = {'3'};
% A={'1';'2';'3';'4';'5';'6';'7';'8';'9';'10';'11';'12';'13';'14';'15';'16';'17';'18';'19';'20';'21';'22'};
for j = 1:22
    i = A{j};
    disp(i)
    bnum = 400;cont=1;
    %output_filename = [outprefix 'Chr' i '.probhap.out'];% Target Haplotype
    %Hapcut_frag1.probhap.corrected.out
    output_filename = [outprefix 'Hapcut_frag' i '.probhap.corrected.out'];
    fid_h=fopen(output_filename,'r');%
    RR=0;
    time=0;
    while(1)
%         if bnum==416 || bnum==513
%             bnum=bnum+1;
%             continue
%         end

        %21
%         if bnum==423 || bnum==513
%             bnum=bnum+1;
%             continue
%         end
        tStart=tic;     %calculate CPU time
        input_filename = [prefix 'Chr' i '_block' num2str(bnum) '.A'];% SNP matrix (0,1,-1) format
        if exist(input_filename,'file')
            %[Matrix,Orginal]=Transform2Geraci(input_filename);
            [Matrix,Orginal]=Transform2Geraci(input_filename);
            [R,C]=size(Matrix);
            %++++++++++++++++++++++++++++++++++++++++++++++++++++
            target='';
            tline=fgetl(fid_h);
            blockend=strfind(tline,'***');
            while isempty(blockend)
                index=strfind(tline,'positions');
                tmp=tline(index+10:length(tline));
                L=str2num(tline(index+10:length(tline)));
                if L==C
                    inthaps=[0,0];
                    r1=1;
                    tline=fgetl(fid_h);
                    blockend=strfind(tline,'***');
                    while isempty(blockend)
                        measures=str2num(tline);
                        if length(measures)==6
                            inthaps(r1,:)=measures(2:3);
                            r1=r1+1;
                        else
                            inthaps(r1,:)=[-1,-1];
                            r1=r1+1;
                        end
                        tline=fgetl(fid_h);
                        blockend=strfind(tline,'***');
                    end
                    numbers=inthaps';
                    %numbers=str2num(tline);%%%%%%%%%%%
                    H1='';H2='';
                    [~,L]=size(numbers);
                    for i2=1:L
                        if numbers(1,i2)==1
                            H1=[H1,'a'];
                        elseif numbers(1,i2)==0
                            H1=[H1,'t'];
                        elseif numbers(1,i2)==-1
                            H1=[H1,'-'];
                        end
                        if numbers(2,i2)==1
                            H2=[H2,'a'];
                        elseif numbers(2,i2)==0
                            H2=[H2,'t'];
                        elseif numbers(2,i2)==-1
                            H2=[H2,'-'];
                        end
                    end
                    target=[H1;H2];
                else
                    'dddddddddddddddddd';
                    tline=fgetl(fid_h);
                    blockend=strfind(tline,'***');
                    while isempty(blockend)
                        tline=fgetl(fid_h);
                        blockend=strfind(tline,'***');
                    end
                end
                
            end
            
            %++++++++++++++++++++++++++++++++++++++++++++++++++++
            if L==C
                    [m0,n0]=size(Matrix);
                    if sizefrag<n0
                        sizefrag=n0
                    end
                    [RR(cont),MECscore(cont)]=NCMHap_realdata(Matrix,target,bnum,Matrix);
                    time(cont)=toc(tStart);% calculate CPU time
                        path='F:\znu\thesis\thesis\code\NCMDiploid\NCMforGithub\1000genome_RR\';
                        fidrr=fopen([path,'RR_chr__' i '.txt'],'w+');%closed
                        for j2=1:cont
                                fprintf(fidrr,'%d\t%1.4f\t%1.4f\t%1.4f\n',j2,RR(j2),MECscore(j2),time(j2));
                        end
                        average_rr=mean(RR);
                        fprintf(fidrr,'meanRR=%1.4f\t',average_rr);
                        average_MEC=mean(MECscore);
                        fprintf(fidrr,'meanMEC=%1.4f\t',average_MEC);
                        average_time=mean(time);
                        fprintf(fidrr,'time=%1.4f',average_time);
                        fclose(fidrr);
                        cont=cont+1;
                        sizefrag
            end
        else
            break
        end
        bnum = bnum + 1         
    end
    fclose(fid_h);
end
