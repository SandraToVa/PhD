load datos.txt
whos r y 
plot(r,y,'o');

%equa='((0.061/r)*(1+b*r+(-((11*pi)/(12*0.061))*(((11*pi/(12*0.061))-1)/(a^2+a*b+(0.187/0.061))))*(r^2))/(1+a*r+(-((11*pi/(12*0.061))-1)/(a^2+a*b+(0.187/0.061)))*(r^2)))+0.187*r+0.559+0.061*(a-b)'
%startPoints = [0.01 0.01]

myfittype = fittype('((0.061/r)*(1+b*r+(-((11*pi)/(12*0.061))*(a^2+a*b+0.187/0.061)/((11*pi)/(12*0.061)-1))*(r^2))/(1+a*r+(-(a^2+a*b+0.187/0.061)/((11*pi)/(12*0.061)-1))*(r^2)))+0.187*r+0.559+0.061*(a-b)','dependent',{'y'},'independent',{'r'},'coefficients',{'a','b'});

%V_p=(0.061/r)*((1+b*r+()*r^2)/(1+a*r+()*r^2))+0.187*r+()
%e=0.559+0.061*(a-b)
%a_2=-(a^2+a*b+0.187/0.061)/((11*pi)/(12*0.061)-1)
%b_2=-((11*pi)/(12*0.061))*(a^2+a*b+0.187/0.061)/((11*pi)/(12*0.061)-1);

%=-46.20097/(a^2+a*b+(0.187/0.061))
%=-2181.548472/(a^2+a*b+(0.187/0.061))

myfit = fit(r',y',myfittype) 
%'Start',startPoints

plot(myfit,r,y);

ci = confint(myfit)


