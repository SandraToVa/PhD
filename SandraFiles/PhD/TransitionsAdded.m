
%IMPORTANT: This code is updated (only the quarkonium part) with the 
% inclusion of the new operator


%As I don't know the best way of making this type of code. This is simply a
%document that acts as a gide of what I have to use inside the loop of
%ComputationExpVal.m

% Uses the functions from FormFactor_ItoF.m and gives the final function
% to compute the transition
%Example: in l=1->l'=1 there are multiple transitions from m->m' that apear
%in I_thetaFuntions. Here we add and obtain the final transitions l->l' in
%terms of funcions m->m'

%This functions will be called in ComputationExpVal in order to directly
%compute the I_ItoF (\, c, s) for each M

%Create a function to call all the functions in here:
function res = TransitionsAdded(x)
    disp(['Input value: ', x]);

    %%%%%%%%%% QQ %%%%%%%%%%%%
    %s->s
    if strcmp(x, 'QQStoS')
        res = @StoStrans;
    %p->p
    elseif strcmp(x, 'QQPtoP')
        res = @PtoPtrans;
    %d->s
    elseif strcmp(x, 'QQDtoS')
        res = @DtoStrans;
    %d->d
    elseif strcmp(x, 'QQDtoD')
        res = @DtoDtrans;
   
    %%%%%%%%%%%% HQ %%%%%%%%%%%%%%
    % p0 -> s
    elseif strcmp(x, 'HQp0tS')
        res = @P0toStrans;
    % p0 -> d
    elseif strcmp(x, 'HQp0tD')
        res = @P0toDtrans;

    % p1 -> s
    elseif strcmp(x, 'HQp1tS') 
        res = @P1toStrans;
    % p1 -> d
    elseif strcmp(x, 'HQp1tD')
        res = @P1toDtrans;

    % (s/d)1 -> p
    elseif strcmp(x, 'HQsdtP') 
        res = @SDtoPtrans;

    % (p/f)2 -> s
    elseif strcmp(x, 'HQpftS')
        res = @PFtoStrans;
    % (p/f)2 -> d
    elseif strcmp(x, 'HQpftD')
        res = @PFtoDtrans;
    end
    disp(['Output function handle: ', func2str(res)]);
end

%First all the functions in Quarkonium -> Quarkonium +2Pions transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% s->s
function StoS=StoStrans
%només 1 contribució Ji=0,m=0->Jf=0,m'=0
%al spin average contribueixen els termes I^2, I_c^2 i I*I_c

%canviar
%estats inciials i finals
Ni=n;
Nf=n';

%m=0
transition0=FormFactor_ItoF('QQS0toS0_F/');
transitionC0=FormFactor_ItoF('QQS0toS0_Fc');

ExpValFunc=ExpValFunctions('QQ');

I_if0(n)=ExpValFunc(Ni,Nf,M(n),transition0,0,0);
I_if_c0(n)=ExpValFunc(Ni,Nf,M(n),transitionC0,0,0);

%el spin average = com nomes hi ha 1 m possible per a la transicio
%directament per a cada

I_if_square(n)=I_if0(n)^2;
I_if_c_square(n)=I_if_c0(n)^2;
I_if_0c_square(n)=I_if0(n)*I_if_c0(n);

end

% p->p
function PtoP=PtoPtrans
% 3 contributions: m=0->m'=0, m=1->m'=1 (estes dos son la mateixa)
% m=-1->m'=-1 for the I, Ic form factors. We will ned a spin average
% (sumatori sobre estats finals i average sobre inicals)

% 1 contribution: m=+1->m'=-1 for the Ix and the same m=-1->m'=+1 for the
% Is case

%canviar
%estats inciials i finals
Ni=n;
Nf=n';

%m=0
transition0=FormFactor_ItoF('QQP0toP0_F/');
transitionC0=FormFactor_ItoF('QQP0toP0_Fc');
%m=1
transition1=FormFactor_ItoF('QQP1toP1_F/');
transitionC1=FormFactor_ItoF('QQP1toP1_Fc');
transitionS1=FormFactor_TtoF('QQP1toP1_Fs');

ExpValFunc=ExpValFunctions('QQ');

I_if0(n)=ExpValFunc(Ni,Nf,M(n),transition0,1,1);
I_if_c0(n)=ExpValFunc(Ni,Nf,M(n),transitionC0,1,1);
I_if1(n)=ExpValFunc(Ni,Nf,M(n),transition1,1,1);
I_if_c1(n)=ExpValFunc(Ni,Nf,M(n),transitionC1,1,1);
I_if_s1(n)=ExpValFunc(Ni,Nf,M(n),transitionS1,1,1);

%el spin average

I_if_square(n)=(I_if0(n)^2 + 2*I_if1(n)^2)/3; % *2 because trans from +-1 to +-1
I_if_c_square(n)=(I_if_c0(n)^2 + 2*I_if_c1(n)^2)/3;
I_if_0c_square(n)=(I_if0(n)*I_if_c0(n) + 2*I_if1(n)*I_if_c1(n))/3;
I_if_s_square(n)=(2*I_if_s1(n)^2)/3;  % *2 because trans from +-1 to -+1

end

% d->s
function DtoS=DtoStrans
% 1 contribution: m=0->m'=0 for the I, Ic form factors

% 1 contribution: m=+2->m'=0 for the Ix and the same m=-2->m'=0 for the
% Is case

%canviar
%estats inciials i finals
Ni=n;
Nf=n';

%m=0
transition0=FormFactor_ItoF('QQD0toS0_F/');
transitionC0=FormFactor_ItoF('QQD0toS0_Fc');
%m=2
transitionS2=FormFactor_ItoF('QQD2toS0_Fs');

ExpValFunc=ExpValFunctions('QQ');

I_if0(n)=ExpValFunc(Ni,Nf,M(n),transition0,2,0);
I_if_c0(n)=ExpValFunc(Ni,Nf,M(n),transitionC0,2,0);
I_if_s2(n)=ExpValFunc(Ni,Nf,M(n),transitionS2,2,0);

%el spin average = com nomes hi ha 1 m possible per a la transicio
%directament per a cada

I_if_square(n)=(I_if0(n)^2)/5;
I_if_c_square(n)=(I_if_c0(n)^2)/5;
I_if_0c_square(n)=(I_if0(n)*I_if_c0(n))/5;
I_if_s_square(n)=(I_if_s2(n)^2)/5;


end

% d->d
function DtoD=DtoDtrans
% 5 contributions: m=0->m'=0, m=1->m'=1 (estes dos son la mateixa)
% m=-1->m'=-1, m=+-2->m'=+-2 for the I, Ic form factors. We will ned a spin average
% (sumatori sobre estats finals i average sobre inicals)

%canviar
%estats inciials i finals
Ni=n;
Nf=n';

%m=0
transition0=FormFactor_ItoF('QQD0toD0_F/');
transitionC0=FormFactor_ItoF('QQD0toD0_Fc');
%m=1
transition1=FormFactor_ItoF('QQD1toD1_F/');
transitionC1=FormFactor_ItoF('QQD1toD1_Fc');
%m=1
transition2=FormFactor_ItoF('QQD2toD2_F/');
transitionC2=FormFactor_ItoF('QQD2toD2_Fc');

ExpValFunc=ExpValFunctions('QQ');

I_if0(n)=ExpValFunc(Ni,Nf,M(n),transition0,2,2);
I_if_c0(n)=ExpValFunc(Ni,Nf,M(n),transitionC0,2,2);
I_if1(n)=ExpValFunc(Ni,Nf,M(n),transition1,2,2);
I_if_c1(n)=ExpValFunc(Ni,Nf,M(n),transitionC1,2,2);
I_if2(n)=ExpValFunc(Ni,Nf,M(n),transition2,2,2);
I_if_c2(n)=ExpValFunc(Ni,Nf,M(n),transitionC2,2,2);

%el spin average

I_if_square(n)=(I_if0(n)^2 + 2*I_if1(n)^2 + 2*I_if2(n)^2)/5; % *2 because trans from +-1 to +-1 and from +-2 to +-2
I_if_c_square(n)=(I_if_c0(n)^2 + 2*I_if_c1(n)^2 + 2*I_if_c2(n)^2)/5;
I_if_0c_square(n)=(I_if0(n)*I_if_c0(n) + 2*I_if1(n)*I_if_c1(n) + 2*I_if2(n)*I_if_c2(n))/5;

end


%First all the functions in Quarkonium -> Quarkonium +2Pions transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Here we take into acount the partial- and partial+ parts of the lagrangian
%density
%There will be cases where we sum over before square because it is the same
%transition

% p0->s
function P0toS=P0toStrans
%Only 1 M->m'
%same contribution from partial- and partial+: multiplied x2
%sum over final and avergae over intial j=0

transition0=I_thetaFunctions('HQp0toS0');

ExpValFunc=ExpValFunctions('HQ');

Itheta0(n)=ExpValFunc(N,N',M(n),transition0,0,0,false);

Itheta(n)=(2*Itheta0(n))*conj(2*Itheta0(n));

end

% p0->d
function P0toD=P0toDtrans
%Only 1 M->m'
%same contribution from partial- and partial+: multiplied x2
%sum over final and avergae over intial j=0

transition0=I_thetaFunctions('HQp0toD0');

ExpValFunc=ExpValFunctions('HQ');

Itheta0(n)=ExpValFunc(N,N',M(n),transition0,0,2,false);

Itheta(n)=(2*Itheta0(n))*conj(2*Itheta0(n));

end

% p1->s
function P1toS=P1toStrans
%Only 1 M->m'
%same contribution from partial- and partial+: diferent sign so =0
%sum over final and avergae over intial j=0

Itheta(n)=0;

end

% p1->d
function P1toD=P1toDtrans
%3 M->m': the one with M=0->m'=0 gives =0 but the others do not
%%so we have M
%same contribution from partial- and partial+: diferent sign so =0
%sum over final and avergae over intial j=0

%M=0->m'=0 when adding partial- and partial+ =0
%M=-1->m'=-1 same but with diferent sign as M=1->m'=1 so they are summed
%over squared

transition1=I_thetaFunctions('HQp1toD1');

ExpValFunc=ExpValFunctions('HQ');

Itheta1(n)=ExpValFunc(N,N',M(n),transition1,1,2,false);

%The N is not n in np_J: remember the ordering of states for (s/d)1 and p1

Itheta(n)=(2*(Itheta1(n)*conj(Itheta1(n))))/3;

end

% (s/d)1->p
function SDtoP=SDtoPtrans
%3 M->m': partial- and partial+ are the same so nothing is 0

%M=0->m'=0 when adding partial- and partial+ = 0
%For partial-: M=-1->m'=-1 is the same x(-1) that for partial+: is M=1->m'=1
%For partial-: M=1->m'=1 is the same x(-1) that for partial+: is M=-1->m'=-1 

%Therefore the transitions 1->1 and -1->-1 give the same result

transition0=I_thetaFunctions('HQsdtoP0');
transition1=I_thetaFunctions('HQsdtoP1');
transitionN=I_thetaFunctions('HQsdtoPn'); %n of negative

ExpValFunc=ExpValFunctions('HQ');

Itheta0(n)=ExpValFunc(N,N',M(n),transition0,1,1,true);
Itheta1(n)=ExpValFunc(N,N',M(n),transition1,1,1,true);
IthetaN(n)=ExpValFunc(N,N',M(n),transitionN,1,1,true);

%The N is not n in n(s/d)_J: remember the ordering of states for (s/d)1 and p1

Itheta(n)=( 2*(Itheta1(n)-IthetaN(n))*conj(Itheta1(n)-IthetaN(n)) )/3;

end

% (p/f)2->s
function PFtoS=PFtoStrans
%1 M->m': partial- and partial+ are the same so nothing is 0

transition0=I_thetaFunctions('HQpftoS0');

ExpValFunc=ExpValFunctions('HQ');

Itheta0(n)=ExpValFunc(N,N',M(n),transition0,2,0,true);

%The N is not n in n(p/f)_J: remember the ordering of states for (p/f)2 and
%d2

Itheta(n)=((2*Itheta0(n))*conj(2*Itheta0(n)))/5;

end

% (p/f)2->s
function PFtoD=PFtoDtrans
%5 M->m': partial- and partial+ are the same so nothing is 0

%M=0->m'=0 when adding partial- and partial+ its the double (all before
%squaring)
%For partial-: M=-1,-2->m'=-1,-2 is the same that for partial+: is M=1,2->m'=1,2 
%For partial-: M=1,2->m'=1,2 is the same that for partial+: is M=-1,-2->m'=-1,-2 

%Therefore the transitions 1->1 and -1->-1; 2->2 and -2->-2 give the same result

transition0=I_thetaFunctions('HQpftoD0');
transition1=I_thetaFunctions('HQpftoD1');
transitionN=I_thetaFunctions('HQpftoDn'); %n of negative
transition2=I_thetaFunctions('HQpftoD2');
transitionM=I_thetaFunctions('HQpftoDm'); %m of doble negative

ExpValFunc=ExpValFunctions('HQ');

Itheta0(n)=ExpValFunc(N,N',M(n),transition0,2,2,true);
Itheta1(n)=ExpValFunc(N,N',M(n),transition1,2,2,true);
IthetaN(n)=ExpValFunc(N,N',M(n),transitionN,2,2,true);
Itheta2(n)=ExpValFunc(N,N',M(n),transition2,2,2,true);
IthetaM(n)=ExpValFunc(N,N',M(n),transitionM,2,2,true);

%The N is not n in n(p/f)_J: remember the ordering of states for (p/f)2 and
%d2

Itheta(n)=(2*(Itheta1(n)+IthetaN(n))*conj(Itheta1(n)+IthetaN(n)) + 2*(Itheta2(n)+IthetaM(n))*conj(Itheta2(n)+IthetaM(n)) + (2*Itheta0(n))*conj(2*Itheta0(n)))/5;

end
