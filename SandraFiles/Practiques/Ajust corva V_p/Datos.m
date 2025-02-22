x = .02:.00001:14;
% En r=.02:.00001:14 dona un rang més proxim a la zona de interès
y=(0.061./x).*(1+0.06964*x-1.45934*x.^2)./(1-0.06733*x+0.01433*x.^2)+0.187.*x+0.551;
A = [x;y];
    
fileID = fopen('datos.txt','w');
%fprintf(fileID,'%6s %12s\n','x','Ruben');
fprintf(fileID,'%12.8f %12.8f\n',A);
fclose(fileID);


type datos.txt
 
