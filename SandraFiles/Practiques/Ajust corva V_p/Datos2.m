r = .02:.00001:14;
f=(0.061./r).*((1-((0.187./0.061)./(1-(12.*0.061)./(11*pi))).*r.^2)./(1-((0.187./0.061)./((11.*pi)./(12.*0.061)-1)).*r.^2))+0.187.*r+0.559;
B = [r;f];
    
fileID = fopen('datos2.txt','w');
%fprintf(fileID,'%6s %12s\n','x','Ruben');
fprintf(fileID,'%12.8f %12.8f\n',B);
fclose(fileID);


type datos2.txt