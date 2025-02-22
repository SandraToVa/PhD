load datos.txt
plot(x,y,'o');

startpoints=[0];
ft = fittype('(0.061/x)*((1+b*x-((0.187/0.061)/(1-(12*0.061)/(11*pi)))*x^2)/(1-((0.187/0.061)/((11*pi)/(12*0.061)-1))*x^2))+0.187*x+(0.559-0.061*b)');
[curve,gof] = fit(x',y',ft, 'Start',startpoints)

hold on
plot(curve,'m')
legend('Data','Ajust')
hold off


%V_p=(0.061/x)*((1+b*x+()*x^2)/(1+()*x^2))+0.187*x+()
%e=0.559-0.061*b
%a_2=-(0.187/0.061)/((11*pi)/(12*0.061)-1)
%b_2=-(0.187/0.061)/(1-(12*0.061)/(11*pi));

license