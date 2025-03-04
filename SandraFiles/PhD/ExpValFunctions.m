%Different functions for QQ and HQ and HH because the operators i diferent
%For Quarkonium-Quarkonium the operators used is one because the wave
%functions are all 1 row vectors
%For Hybrid-Quarkonium the operators are 2 because the hybrids have 2 rows
%and for each row the operator is diferent so we need to import 2 operators


function res = ExpValFunctions(x)
    if strcmp(x, 'QQ')
        res = @QQComputeExpValM;
    elseif strcmp(x, 'HQ')
        res = @HQComputeExpValM;
    end
end



function EV=QQComputeExpValM(i,f,M,trans,Ji,Jf)
%With a initial al final states we compute expectations values
%The operators of the expectation depend on M (dipion mass)
%We want to find the relationship between the result and the M used
%trans=transition function StoS, PtoP, etc

%initail
setspin(0);
sin=spin;
sfin=spin;

%Define the funtions to use
QuarkoniumS0Ji=str2func(['QuarkoniumS0J' num2str(Ji)]);
QuarkoniumS0Jf=str2func(['QuarkoniumS0J' num2str(Jf)]);

[Einicial,~,xi]=QuarkoniumS0Ji(m_q,sin);

Ei=Einicial(i);

%final
[Efinal,~,xf]=QuarkoniumS0Jf(m_q,sfin);

Ef=Efinal(f);

E1=abs(Ef-Ei);
if M>E1
    EV=0;
    disp(['M>E1 cannot be', num2str(E1)])

else
%Operator that we want to find the expectated value
%The x must be the same. It is a problem of tolerance if not
    if length(xi)~=length(xf)
     disp('Change tolerance:');
     disp(length(xi)-length(xf));
    else
     x=xi; %The same for xi and xf same x
     Op=trans(x,Ei,Ef,M);
     Op=diag(diag(Op));

     %PROBLEM WITH THE FUNCTION, x for hybrid and x for quarkonium have diferent
     %dimension. Due to ConstrucMesh?????
     V=ExpectedValue(i,f,Op,QuarkoniumS0Ji,QuarkoniumS0Jf,m_q,sin,sfin);

     %Si sabem les funcions d'ona pq ens les donen i les energies utilitzem
     %la seguent funcio
     %V=ExpectedValueKnownWF(Wi,Wf,x,Op,s)

     %Only one value is diferent from zero
     EV=V;
     %[row,col]=find(V);
     %EV=V(row,col);
    end
end

end

function EV=HQComputeExpValM(i,f,M,trans,Ji,Jf,double)
%With a initial al final states we compute expectations values
%The operators of the expectation depend on M (dipion mass)
%We want to find the relationship between the result and the M used
%trans=transition function from I_thetaFunctions respository
%double= if the initial state (the hybrid) is a state composed by two
%(s/d)1 or (p/f)2 double is true and the operator has 2 components and the
%w.f. is not only one like but two

%Define the funtions to use
QuarkoniumS0Ji=str2func(['QuarkoniumS0J' num2str(Ji)]);
QuarkoniumS0Jf=str2func(['QuarkoniumS0J' num2str(Jf)]);

%initail
setspin(1);
sin=spin;

[Einicial,~,xi]=QuarkoniumS0Ji(m_q,sin);

Ei=Einicial(i);

%final
setspin(0)
sfin=spin;

[Efinal,~,xf]=QuarkoniumS0Jf(m_q,sfin);

Ef=Efinal(f);

E1=abs(Ef-Ei);
if M>E1
    EV=0;
    disp(['M>E1 cannot be', num2str(E1)])

else
%Operator that we want to find the expectated value
%The x must be the same. It is a problem of tolerance if not
    if length(xi)~=length(xf)
     disp('Change tolerance:');
     disp(length(xi)-length(xf));
    else
    %We must direnciate for hybrids with 2component wf and hybrids with
    %1component wf
    %create a variable: "double". If true, it has double wf. If false it has
    %not
     x=xi; %The same for xi and xf same x

     if double
         %Double hybrids have 3 rows and the 1st and 2nd are treated
         %searately and the 3rd is zeros

         %Important: The 1st row is the one that will multiply the 1st row
         %of the wf. They are s from (s/d)1 and p from (p/f)2
         Op_d=trans(x,Ei,Ef,M);  %double
         Op_d1=Op_d(:,:,1);
         Op_d2=Op_d(:,:,2);
         [numR,numC] = size(Op_d1);

         Op(:,:,1)=diag(diag(Op_d1));
         Op(:,:,2)=diag(diag(Op_d2));
         Op(:,:,3)=zeros(numR,numC);


     else
         %Non double hybrids have 3 rows but only the third/first has
         %information (or one or the other)
         Op_s=trans(x,Ei,Ef,M);   %single
         [numR,numC] = size(Op_s);
         Op_s=diag(diag(Op_s));
         Op(:,:,1)=Op_s;
         Op(:,:,2)=zeros(numR,numC);
         Op(:,:,3)=Op_s;
     end

     %PROBLEM WITH THE FUNCTION, x for hybrid and x for quarkonium have diferent
     %dimension. Due to ConstrucMesh?????
     V=ExpectedValue(i,f,Op,QuarkoniumS0Ji,QuarkoniumS0Jf,m_q,sin,sfin);

     %SHA DE CANVIAR PER A FER LA SUMA SI DOUBLE I NO SI SINGLE!! o meibi no
     %fa falta
     EV=sum(V);

    end
end

end


%VALOR DE LA MASSA
%m=1.4702; charm
%m=4.8802; bottom
function x0=m_q
global v0
x0=v0;
end

function setm_q(val0)
global v0
v0=val0;
end

% VALOR DEL SPIN
function x7=spin
global v7
x7=v7;
end 

function setspin(val7)
global v7
v7 = val7;
end

% VALOR DEL r0
function x3=r0
global v3
x3=v3;
end 

function setr0(val3)
global v3
v3 = val3;
end

