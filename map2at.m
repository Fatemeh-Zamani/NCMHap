function str=map2at(inStr)
inStr=char(inStr);
L=length(inStr);
AT='';
for i=1:L
    if inStr(i)=='1'
        AT=[AT,'a'];        
    elseif inStr(i)=='0'
        AT=[AT,'t'];        
    else
        AT=[AT,'-'];        
    end
end
str=AT;
end