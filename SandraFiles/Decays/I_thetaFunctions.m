

%CODE NOT VALID. The new version is FormFactors_ItoF




% UPDATE: En lloc de Itheta son Ii->f. And the sign is diferent! Instead of
% - with e^{-iphi} and + with e^{+iphi} it is - for both
% Només ho he canviat per a (s/d)1

%All the functions for computing I_\theta 
%These functions are the angular dependent part (theta and fi) integrated
%with the corresponendt spherical harmonics of the states we want to
%compute the decay with of.
%Actualitzat en la integracio de 0 a pi 

%Create a function to call all the functions in here:
function res = I_thetaFunctions(x)
    disp(['Input value: ', x]);
    %%%%%%%%%% QQ %%%%%%%%%%%%
    if strcmp(x, 'QQS0toS0')
        res = @StoStrans;
    %elseif strcmp(x, 'QQP0toS0') %No possible =0
    %    res = @P0toStrans;
    elseif strcmp(x, 'QQP0toP0')
        res = @P0toP0trans;
    elseif strcmp(x, 'QQP1toP1')
        res = @P1toP1trans;
    elseif strcmp(x, 'QQD0toS0')
        res = @DtoStrans;
    %elseif strcmp(x, 'QQD0toP0') %No possible =0
    %    res = @D0toP0trans;
    %elseif strcmp(x, 'QQD1toP1') %No possible =0
    %    res = @D1toP1trans;
    elseif strcmp(x, 'QQD0toD0')
        res = @D0toD0trans;
    elseif strcmp(x, 'QQD1toD1')
        res = @D1toD1trans;
    elseif strcmp(x, 'QQD2toD2')
        res = @D2toD2trans;
   
    %%%%%%%%%%%% HQ %%%%%%%%%%%%%%
    % p0 trans
    elseif strcmp(x, 'HQp0toS0')
        res = @Hp0toQs0trans;
    elseif strcmp(x, 'HQp0toD0')
        res = @Hp0toQd0trans;
    % p1 trans
    elseif strcmp(x, 'HQp1toS0') 
        res = @Hp1toQs0trans;
    elseif strcmp(x, 'HQp1toD0')
        res = @Hp1toQd0trans;
    elseif strcmp(x, 'HQp1toD1')
        res = @Hp1toQd1trans;
    % (s/d)1 trans
    elseif strcmp(x, 'HQsdtoPn') %m'=-1 = n
        res = @HsdtoQpntrans;
    elseif strcmp(x, 'HQsdtoP0')
        res = @HsdtoQp0trans;
    elseif strcmp(x, 'HQsdtoP1')
        res = @HsdtoQp1trans;
    % (p/f)2 trans
    elseif strcmp(x, 'HQpftoS0')
        res = @HpftoQs0trans;
    elseif strcmp(x, 'HQpftoDm') %m'=-2 = m
        res = @HpftoQdmtrans;
    elseif strcmp(x, 'HQpftoDn') %m'=-1 = n
        res = @HpftoQdntrans;
    elseif strcmp(x, 'HQpftoD0')
        res = @HpftoQd0trans;
    elseif strcmp(x, 'HQpftoD1')
        res = @HpftoQd1trans;
    elseif strcmp(x, 'HQpftoD2')
        res = @HpftoQd2trans;
    end
    disp(['Output function handle: ', func2str(res)]);
end



%First all the functions in Quarkonium -> Quarkonium +2Pions transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%With the functions below we tant to create a matrix "operator"
%This matrix is diagonal, however I sum and substract elements to each
%component of the matrix resuylting in a non diagonal matrix.
%In the function ComputeExpVal I only use the matrix elements of these
%matrix and hence it does not matter that the matrix from the functions are
%not diagonal
%But IMPORTANT to remember!!!

%For s->s transitions
function Fss=StoStrans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%Fss=r/(4*pi)-(r^3*Sq^2)/(288*pi);
%Hole function
Fss=(1/(2*pi*Sq))*sinint((Sq/2).*r);
end

%For p->s transitions (l=1,m=0 -> l=0)
function Fps=P0toStrans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Fps=0.*r;
end

%For p->p transitions (l=1,m=0 -> l=1,m=0)
function Fpp00=P0toP0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
Fpp00=(3/pi)*(2*sin((Sq/2).*r)-(Sq.*r).*cos((Sq/2).*r))/(Sq^3.*r^2);
end

%For p->p transitions (l=1,m=+-1 -> l=1,m=+-1)
function Fpp11=P1toP1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Fpp11=(3/(4*pi*Sq))*(-(sin((Sq/2).*r)-((Sq/2).*r).*cos((Sq/2).*r))/((Q/4).*r^2)+sinint((Sq/2).*r));
end

%For d->s transitions (l=2,m=0 -> l=0)
function Fds=DtoStrans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
%Aproximation
%Fds=-(r^3*Sq^2)/(144*(sqrt(5)*pi));
%Hole function
Fds=-(sqrt(5)/(8*pi))*(2*sinint((Sq/2).*r)/Sq+3*(4*Sq*cos((Sq/2).*r).*r - 8*sin((Sq/2).*r))/(Sq^3.*r^2));
end

%For d->p transitions (l=2,m=0 -> l=1,m=0)
function Fdp00=D0toP0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Fdp00=0.*r;
end

%For d->p transitions (l=2,m=+-1 -> l=1,m=+-1)
function Fdp11=D1toP1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Fdp11=0.*r;
end

%For d->d transitions (l=2,m=0 -> l=2,m=0)
function Fdd00=D0toD0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%Aprox
%Fdd00=r/(4*pi)-(11*r^3*Q)/(2016*pi);
%Hole
Fdd00=(5/(8*pi*Sq))*(sinint((Sq/2).*r)+6*(-Sq.*r*(Q.*r^2-72).*cos((Sq/2).*r)+2*(7*Q.*r^2-72).*sin((Sq/2).*r))/(Q^2.*r^4));
end

%For d->d transitions (l=2,m=+-1 -> l=2,m=+-1)
function Fdd11=D1toD1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%Aprox
%Fdd11=r/(4*pi)-(r^3*Q)/(224*pi);
%Hole
Fdd11=(-30/pi)*((6*Sq.*r.*cos((Sq/2).*r)+(Q.*r^2 - 12).* sin((Sq/2).*r))/(Sq^5.*r^4));
end

%For d->d transitions (l=2,m=+-2 -> l=2,m=+-2)
function Fdd22=D2toD2trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%Aprox
%Fdd22=r/(4*pi)-(r^3*Q)/(672*pi);
%Hole
Fdd22=(15/(16*pi*Sq))*(sinint((Sq/2).*r)+(2*Sq.*r*(Q.*r^2+24).*cos((Sq/2).*r)+4*(Q.*r^2-24).*sin((Sq/2).*r))/(Q^2.*r^4));
end

%Now we start with the functions Hybrid -> Quarkonium +2Pions transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We must remember that all these are multiplied by an imaginary (*1i)
%But as the part we are interested in is the product with the complex
%conjugate the imaginary part plays no rol.

% IMPORTANT - WRONG!!
%The only functions actually used are the (s/d)1 functions. For them the
%behaviour was wrong unless using the tyalor development at Sq*r small.
%This is missing in the other functions. If we want to use them we should
%implement the taylor development

%%%%%%%%% p0  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Transition between Hp0->Q(l=0,m=0): HQp0toS0
%Sandwitch directe
function Hp0s0=Hp0toQs0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%Green integral
A= 1i / (8 * pi);
B=diag(diag((1 ./ (Q * r * pi))));
C1=4 * pi * sin((Sq / 2) * r);
C3=(pi - Sq * r) * (pi + Sq * r) ;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C3*C4;
Hp0s0=A*B*C;
end

%Transition between Hp0->Q(l=2,m=0): HQp0toD0
%Sandwitch directe
function Hp0d0=Hp0toQd0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%Blue integral
A= 1i * sqrt(5) / (16 * pi);
B=diag(diag((1 ./ (Q^2 * r^3 * pi))));
C1=48 * pi * Sq * r * cos((Sq / 2) * r);
C2=4 * pi * (-24 + 3*pi^2 - Q*r^2) * sin((Sq / 2) * r);
C3=3*pi^4 - 4*pi^2*Q*r^2 + Q^2*r^4;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hp0d0=A*B*C;
end


%%%%%%%%%%% p1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Transition between Hp1->Q(l=0,m=0): HQp1toS0
%Only contribution for Q(l=0)
%Sandwitch directe
function Hp1s0=Hp1toQs0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
A= 1i * sqrt(6) / (16 * pi);
B=diag(diag((1 ./ (Q * r * pi))));
C1=4 * pi * sin((Sq / 2) * r);
C3=(pi - Sq * r) * (pi + Sq * r) ;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C3*C4;
Hp1s0=A*B*C;
end

%Transition between Hp1->Q(l=2,m=0): HQp1toD0
%Two contribution for Q(l=2)
%Sandwitch directe
function Hp1d0=Hp1toQd0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
A= 1i * sqrt(3) / (32 * pi);
B=diag(diag((1 ./ (Q^2 * r^3 * pi))));
C1=48 * pi * Sq * r * cos((Sq / 2) * r);
C2=4 * pi * (-24 + 3*pi^2 - Q*r^2) * sin((Sq / 2) * r);
C3=3*pi^4 - 4*pi^2*Q*r^2 + Q^2*r^4;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hp1d0=A*B*C;
end

%Transition between Hp1->Q(l=2,m=-1): HQp1toD0
%Two contribution for Q(l=2)
%Sandwitch directe
function Hp1d1=Hp1toQd1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
A= -1i * 3 * sqrt(10) / (16 * pi);
B=diag(diag((1 ./ (Q^2 * r^3))));
C1=16 * Sq * r * cos((Sq / 2) * r);
C2=4 * (-8 + pi^2) * sin((Sq / 2) * r);
C3=pi * (pi - Sq*r) * (pi + Sq*r);
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hp1d1=A*B*C;
end


%%%%%%%%%% (s/d)1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Transition between Hsd1->Q(l=1,m=-1): HQsdtoP-
%3 contribution for Q(l=1) this is one of them
%Sandwitch: 1st row L=0 (s) and 2nd row L=2 (d) part of the hyrbid
function Hsdpn=HsdtoQpntrans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);

% Loop through each element of the diagonal matrix r
[rows, cols] = size(r);
Hsdpn(:,:,:)=zeros(rows,cols,2);
for i = 1:rows
    % Get the diagonal element
    r_elem = r(i, i);
    
    % Compute the product of Sq and the diagonal element
    product = Sq * r_elem;
    
% L=0 (part s del hybrid)
    A= -1i * sqrt(3) / (8 * pi); %-1i *
    % Apply the appropriate function based on the product value
    if true %product < 1/10
        Hsdpn(i,i,1)= A * ( ( (4*r_elem^2*Sq) / (3*pi^2) ) + ( (8 - pi^2)*r_elem^4*Sq^3 / (30*pi^4) ) );
    else
        B=1 / (Q * r_elem * pi);
        C1=4 * pi * sin((Sq / 2) * r_elem);
        C3=(pi - Sq * r_elem) * (pi + Sq * r_elem) ;
        C4=sinint((pi - Sq * r_elem) / 2) - sinint((pi + Sq * r_elem) / 2);
        C=C1+C3*C4;
        Hsdpn(i,i,1)=A*B*C;
    end
    
%L=2 (d part of the hybrid)
     A= -1i * sqrt(6) / (32 * pi); %-1i *
    % Apply the appropriate function based on the product value
    if true %product < 1/10
        Hsdpn(i,i,2)= A * ( -(8*r_elem^2*Sq)/(15*pi^2) + (8-pi^2)*r_elem^4*Sq^3/(105*pi^4) );
    else
        B=1 / (Q^2 * r_elem^3 * pi);
        C1=48 * pi * Sq * r_elem * cos((Sq / 2) * r_elem);
        C2=4 * pi * (-24 + 3*pi^2 - Q*r_elem^2) * sin((Sq / 2) * r_elem);
        C3=(3*pi^4 - 4*pi^2*Q*r_elem^2 + Q^2*r_elem^4);
        C4=sinint((pi - Sq * r_elem) / 2) - sinint((pi + Sq * r_elem) / 2);
        C=C1+C2+C3*C4;
        Hsdpn(i,i,2)=A*B*C;
    end
end %for end

end %function end

%Transition between Hsd1->Q(l=1,m=0): HQsdtoP0
%3 contribution for Q(l=1) this is the 2nd one
%Sandwitch: 1st row L=0 (s) and 2nd row L=2 (d) part of the hyrbid
function Hsdp0=HsdtoQp0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion p0invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);

% Loop through each element of the diagonal matrix r
[rows, cols] = size(r);
Hsdp0(:,:,:)=zeros(rows,cols,2);
for i = 1:rows
    % Get the diagonal element
    r_elem = r(i, i);
    
    % Compute the product of Sq and the diagonal element
    product = Sq * r_elem;
    
% L=0 (part s del hybrid)
%this is zero   
    
%L=2 (d part of the hybrid)
     A= 1i * 3 * sqrt(6) / (16 * pi); %1i *
    % Apply the appropriate function based on the product value
    if true %product < 1/10
        Hsdp0(i,i,2)= A * ( (4*r_elem^2*Sq)/(15*pi^2) + (8-pi^2)*r_elem^4*Sq^3/(70*pi^4) );
    else
        B=1 / (Q^2 * r_elem^3);
        C1=16 * Sq * r_elem * cos((Sq / 2) * r_elem);
        C2=4 * (-8 + pi^2) * sin((Sq / 2) * r_elem);
        C3=pi * (pi - Sq*r_elem) * (pi + Sq*r_elem);
        C4=sinint((pi - Sq * r_elem) / 2) - sinint((pi + Sq * r_elem) / 2);
        C=C1+C2+C3*C4;
        Hsdp0(i,i,2)=A*B*C;
    end
end %for end

end %function end


%Transition between Hsd1->Q(l=1,m=1): HQsdtoP1
%3 contribution for Q(l=1) this is the 3rd one
%Sandwitch: 1st row L=0 (s) and 2nd row L=2 (d) part of the hyrbid
function Hsdp1=HsdtoQp1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);

% Loop through each element of the diagonal matrix r
[rows, cols] = size(r);
Hsdp1(:,:,:)=zeros(rows,cols,2);
for i = 1:rows
    % Get the diagonal element
    r_elem = r(i, i);
    
    % Compute the product of Sq and the diagonal element
    product = Sq * r_elem;
    
% L=0 (part s del hybrid)
%this is zero   
    
%L=2 (d part of the hybrid)
     A= 1i * 3 * sqrt(6) / (32 * pi); %1i *
    % Apply the appropriate function based on the product value
    if true %product < 1/10
        Hsdp1(i,i,2)= A *( (16*r_elem^2*Sq)/(15*pi^2) - 2*(-8+pi^2)*r_elem^4*Sq^3/(105*pi^4) );
    else
        B=1 / (Q^2 * r_elem^3 * pi);
        C1=-16 * pi * Sq * r_elem * cos((Sq / 2) * r_elem);
        C2=4 * pi * (8 - pi^2 + Q*r_elem^2) * sin((Sq / 2) * r_elem);
        C3=-(pi^2 - Q*r_elem^2)^2;
        C4=sinint((pi - Sq * r_elem) / 2) - sinint((pi + Sq * r_elem) / 2);
        C=C1+C2+C3*C4; 
        Hsdp1(i,i,2)=A*B*C;
    end
end %for end

end %function end


%%%%%%%%%% (p/f)2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Transition between Hpf2->Q(l=0,m=0): HQpftoS0
%Only contribution for Q(l=0)
%Sandwitch: 1st row L=1 (p) and 2nd row L=3 (f) part of the hyrbid
function Hpfs0=HpftoQs0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=1 (p part of the hybrid)
A= 1i / (8 * pi);
B=diag(diag((1 ./ (Q * r * pi))));
C1=4 * pi * sin((Sq / 2) * r);
C3=(pi - Sq * r) * (pi + Sq * r) ;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C3*C4;
Hpfs0(:,:,1)=A*B*C;

%L=3 (f part of the hybrid)
A=1i * sqrt(3)/ (16 * pi);
B=diag(diag((1 ./ (Q^2 * r^3 * pi))));
C1=80 * pi * Sq * r .* cos((Sq / 2) * r);
C2=4 * pi * (-40 +5*pi^2 - Q*r^2)  .* sin((Sq / 2) * r);
C3=(5*pi^4 - 6*pi^2*Q*r^2 + Q^2*r^4) ;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hpfs0(:,:,2)=A*B*C;
end

%Transition between Hpf2->Q(l=2,m=-2): HQpftoDm
%To Q(l=2) contribute 5 terms (m=-2)
%Sandwitch: 1st row L=1 (p) and 2nd row L=3 (f) part of the hyrbid
function Hpfdm=HpftoQdmtrans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=1 (p part of the hybrid)
A= -1i * 3 * sqrt(10) / (32 * pi);
B=diag(diag((1 ./ (Q^2 * r^3 * pi))));
C1=-16 * pi * Sq * r * cos((Sq / 2) * r);
C2=4 * pi * (8 - pi^2 + Q*r^2) * sin((Sq / 2) * r);
C3=(pi^2 - Q*r^2)^2;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2-C3*C4;
Hpfdm(:,:,1)=A*B*C;
%L=3 (f part of the hybrid)
A= -1i * sqrt(15) / (64 * pi);
B=diag(diag((1 ./ (Q^3 * r^5 * pi))));
C1=16 * pi * Sq * r * (240 - 5*pi^2 + Q*r^2) *cos((Sq / 2) * r);
C2=-4 * pi * (5*(384 - 8*pi^2 + pi^4) - 2*(76 + 3*pi^2)*Q*r^2 + Q^2*r^4) * sin((Sq / 2) * r);
C3=((pi^2 - Q*r^2)^2) * (5*pi^2 - Q*r^2);
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hpfdm(:,:,2)=A*B*C;
end

%Transition between Hpf2->Q(l=2,m=-1): HQpftoDn
%To Q(l=3) contribute 5 terms (m=-1)
%Sandwitch: 1st row L=1 (p) and 2nd row L=3 (f) part of the hyrbid
function Hpfdn=HpftoQdntrans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=1 (p part of the hybrid)
A= -1i * 3 * sqrt(10) / (16 * pi);
B=diag(diag((1 ./ (Q^2 * r^3))));
C1=16 * Sq * r * cos((Sq / 2) * r);
C2=4 * (-8 + pi^2) * sin((Sq / 2) * r);
C3=pi * (pi - Sq*r) * (pi + Sq*r);
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hpfdn(:,:,1)=A*B*C;
%L=3 (f part of the hybrid)
A= -1i * sqrt(15) / (16 * pi);
B=diag(diag((1 ./ (Q^3 * r^5))));
C1=16 * Sq * r * (-240 + 5*pi^2 + 2*Q*r^2) * cos((Sq / 2) * r);
C2=4 * (5*(384 - 8*pi^2 + pi^4) - (176 + 3*pi^2)*Q*r^2) * sin((Sq / 2) * r);
C3=pi * (5*pi^4 - 8*pi^2*Q*r^2 + 3*Q^2*r^4);
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hpfdn(:,:,2)=A*B*C;
end

%Transition between Hpf2->Q(l=2,m=0): HQpftoD0
%To Q(l=2) contribute 5 terms (m=0)
%Sandwitch: 1st row L=1 (p) and 2nd row L=3 (f) part of the hyrbid
function Hpfd0=HpftoQd0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=1 (p part of the hybrid)
A=1i * sqrt(10)/ (32 * pi);
B=diag(diag((1 ./ (Q^2 * r^3 * pi))));
C1=48 * pi * Sq * r * cos((Sq / 2) * r);
C2=4 * pi * (-24 + 3*pi^2 - Q*r^2) * sin((Sq / 2) * r);
C3=3*pi^4 - 4*pi^2*Q*r^2 + Q^2*r^4;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hpfd0(:,:,1)=A*B*C;

%L=3 (f part of the hybrid)
A=1i * sqrt(15)/ (32 * pi);
B=diag(diag((1 ./ (Q^3 * r^5 * pi))));
C1=16 * pi * Sq * r * (-720 + 15*pi^2 + 7*Q*r^2) .* cos((Sq / 2) * r);
C2=4 * pi * (15*(384 - 8*pi^2 + pi^4) - 8*(67 + pi^2)*Q*r^2 + Q^2*r^4)  .* sin((Sq / 2) * r);
C3=(15*pi^6 - 23*pi^4*Q*r^2 + 9*pi^2*Q^2*r^4 - Q^3*r^6) ;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hpfd0(:,:,2)=A*B*C;
end

%Transition between Hpf2->Q(l=2,m=1): HQpftoD1
%To Q(l=2) contribute 5 terms (m=1)
%Sandwitch: 1st row L=1 (p) and 2nd row L=3 (f) part of the hyrbid
function Hpfd1=HpftoQd1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=1 (p part of the hybrid)
Hpfd1(:,:,1)=0 * r;

%L=3 (f part of the hybrid)
A=1i * 5 * sqrt(15)/ (16 * pi);
B=diag(diag((1 ./ (Q^3 * r^5))));
C1=-16 * (-48 + pi^2) * Sq * r * cos((Sq / 2) * r);
C2=4 * (-384 + 8*pi^2 - pi^4 + (32 + pi^2)*Q*r^2) .* sin((Sq / 2) * r);
C3=-pi*(pi^2 - Q*r^2)^2 ;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hpfd1(:,:,2)=A*B*C;
end

%Transition between Hpf2->Q(l=2,m=2): HQpftoD2
%To Q(l=2) contribute 5 terms (m=2)
%Sandwitch: 1st row L=1 (p) and 2nd row L=3 (f) part of the hyrbid
function Hpfd2=HpftoQd2trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=1 (p part of the hybrid)
Hpfd2(:,:,1)=0*r;

%L=3 (f part of the hybrid)
A=1i * 5* sqrt(15)/ (64 * pi);
B=diag(diag((1 ./ (Q^3 * r^5 * pi))));
C1=16 * pi * Sq * r * (-48 + pi^2 - Q*r^2) .* cos((Sq / 2) * r);
C2=4 * pi * (384 - 8*pi^2 + pi^4 - 2*(12 + pi^2)*Q*r^2 + Q^2*r^4)  .* sin((Sq / 2) * r);
C3=(pi^2 - Q*r^2)^3;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hpfd2(:,:,2)=A*B*C;
end
