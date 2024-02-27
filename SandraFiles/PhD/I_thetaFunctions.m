
%All the functions for computing I_\theta
%These functions are the angular dependent part (theta and fi) integrated
%with the corresponendt spherical harmonics of the states we want to
%compute the decay with of.

%Create a function to call all the functions in here:
function res = I_thetaFunctions(x)
    disp(['Input value: ', x]);
    % QQ
    if strcmp(x, 'QQS0toS0')
        res = @StoStrans;
    elseif strcmp(x, 'QQP0toS0')
        res = @P0toStrans;
    elseif strcmp(x, 'QQP0toP0')
        res = @P0toP0trans;
    elseif strcmp(x, 'QQP1toP1')
        res = @P1toP1trans;
    elseif strcmp(x, 'QQD0toS0')
        res = @DtoStrans;
    elseif strcmp(x, 'QQD0toP0')
        res = @D0toP0trans;
    elseif strcmp(x, 'QQD1toP1')
        res = @D1toP1trans;
    elseif strcmp(x, 'QQD0toD0')
        res = @D0toD0trans;
    elseif strcmp(x, 'QQD1toD1')
        res = @D1toD1trans;
    elseif strcmp(x, 'QQD2toD2')
        res = @D2toD2trans;
    % HQ
    elseif strcmp(x, 'HQp0toD2')
        res = @Hp0toQd2trans;
    elseif strcmp(x, 'HQp1toP1')
        res = @Hp1toQp1trans;
    elseif strcmp(x, 'HQp1toD1')
        res = @Hp1toQd1trans;
    elseif strcmp(x, 'HQp1toD2')
        res = @Hp1toQd2trans;
    elseif strcmp(x, 'HQsdtoP1')
        res = @HsdtoQp1trans;
    elseif strcmp(x, 'HQsdtoD1')
        res = @HsdtoQd1trans;
    elseif strcmp(x, 'HQsdtoD2')
        res = @HsdtoQd2trans;
    elseif strcmp(x, 'HQd2toS0')
        res = @Hd2toQs0trans;
    elseif strcmp(x, 'HQd2toP0')
        res = @Hd2toQp0trans;
    elseif strcmp(x, 'HQd2toP1')
        res = @Hd2toQp1trans;
    elseif strcmp(x, 'HQd2toD0')
        res = @Hd2toQd0trans;
    elseif strcmp(x, 'HQd2toD1')
        res = @Hd2toQd1trans;
    elseif strcmp(x, 'HQd2toD2')
        res = @Hd2toQd2trans;
    elseif strcmp(x, 'HQpftoS0')
        res = @HpftoQs0trans;
    elseif strcmp(x, 'HQpftoP0')
        res = @HpftoQp0trans;
    elseif strcmp(x, 'HQpftoP1')
        res = @HpftoQp1trans;
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
Fss=(1/(4*pi*Sq))*sinint((Sq/2).*r);
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
Fps=(sqrt(3)/(4*pi))*(1-cos((Sq/2).*r))/((Q/2).*r);
end

%For p->p transitions (l=1,m=0 -> l=1,m=0)
function Fpp00=P0toP0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
Fpp00=(3/pi)*(sin((Sq/2).*r)-((Sq/2).*r).*cos((Sq/2).*r))/(Sq^3.*r^2);
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
Fpp11=(3/(8*pi*Sq))*(-(sin((Sq/2).*r)-((Sq/2).*r).*cos((Sq/2).*r))/((Q/4).*r^2)+sinint((Sq/2).*r));
end

%For d->s transitions (l=2,m=0 -> l=0)
function Fds=DtoStrans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
r=diag(x);
Fds=-(sqrt(5)/(8*pi))*(sinint((Sq/2).*r)/Sq+6*(Sq*cos((Sq/2).*r).*r-2*sin((Sq/2).*r))/(Sq^3.*r^2));
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
Fdp00=(sqrt(15)/(8*pi))*(-48-2*Q.*r^2+(48-4*Q.*r^2).*cos((Sq/2).*r)+24*Sq.*r.*sin((Sq/2).*r))/(Q^2.*r^3);
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
Fdp11=(3*sqrt(5)/(4*pi))*(8+Q.*r^2-8*cos((Sq/2).*r)-4*Sq.*r.*sin((Sq/2).*r))/(Q^2.*r^3);
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
Fdd00=(5/(16*pi*Sq))*(sinint((Sq/2).*r)+6*(-Sq.*r*(Q.*r^2-72).*cos((Sq/2).*r)+2*(7*Q.*r^2-72).*sin((Sq/2).*r))/(Q^2.*r^4));
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
Fdd11=(-15/pi)*((6*Sq.*r.*cos((Sq/2).*r)+(Q.*r^2 - 12).* sin((Sq/2).*r))/(Sq^5.*r^4));
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
Fdd22=(15/(16*pi*Sq))*(sinint((Sq/2).*r)+2*(Sq.*r*(Q.*r^2+24).*cos((Sq/2).*r)+2*(Q.*r^2-24).*sin((Sq/2).*r))/(Q^2.*r^4));
end

%Now we start with the functions Hybrid -> Quarkonium +2Pions transitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We must remember that all these are multiplied by an imaginary (*1i)
%But as the part we are interested in is the product with the complex
%conjugate the imaginary part plays no rol.

%%%%%%%%% p0
%Transition between Hp0->Q(l=2,m=2): HQp0toD2
%Only transition for p0 states to quarkonium for l=<2
%Sandwitch directe
function Hp0d2=Hp0toQd2trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
A=(sqrt(15) / (32 * pi * sqrt(2)));
B=diag(diag((1 ./ (Q^2 * r^3))));
C1=16 * pi * Sq * r .* cos((Sq / 2) * r);
C2=4 * pi * (-8 + pi^2 - Q * r^2) .* sin((Sq / 2) * r);
C3=(pi^2 - Q * r^2)^2 ;
C4=sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2);
C=C1+C2+C3*C4;
Hp0d2=A*B*C;
%Hp0d2=(sqrt(15) / (32 * pi * sqrt(2))) * (1 ./ (Q^2 * r.^3)) * (16 * pi * Sq * r .* cos((Sq / 2) * r) + 4 * pi * (-8 + pi^2 - Q * r.^2) .* sin((Sq / 2) * r) + (pi^2 - Q * r.^2).^2 .* (sinint((pi - Sq * r) / 2) - sinint((pi + Sq * r) / 2)));
end

%%%%%%%%%%% p1
%Transition between Hp1->Q(l=1,m=1): HQp1toP1
%Only contribution for Q(l=1)
%Sandwitch directe
function Hp1p1=Hp1toQp1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hp1p1=(-3/(16*pi*sqrt(2)))*(1/(Sq^3*r^2))*(-8+8*cos((Sq/2)*r)+4*Sq*r*sin((Sq/2)*r)+(pi^2-Q*r^2)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)-2*sinint(pi/2)));
end

%Transition between Hp1->Q(l=2,m=1): HQp1toD1
%To Q(l=2) contribute 2 terms (m=1)
%Sandwitch directe
function Hp1d1=Hp1toQd1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hp1d1=(-3*sqrt(10)/(32*pi))*(1/(Q^2*r^3))*(16*pi*Sq*r*cos((Sq/2)*r)+4*(-8+pi^2)*sin((Sq/2)*r)+pi*(pi^2-Q*r^2)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
end

%Transition between Hp1->Q(l=2,m=2): HQp1toD2
%To Q(l=2) contribute 2 terms (m=2)
%Sandwitch directe
function Hp1d2=Hp1toQd2trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hp1d2=(3*sqrt(5)/(64*pi^2))*(1/(Q^2*r^3))*(16*pi*Sq*r*cos((Sq/2)*r)+4*pi*(-8+pi^2-Q*r^2)*sin((Sq/2)*r)+(pi^2-Q*r^2)^2*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
end

%%%%%%%%%% (s/d)1
%Transition between Hsd1->Q(l=1,m=1): HQsdtoP1
%Only contribution for Q(l=1)
%Sandwitch: 1st row L=0 (s) and 2nd row L=2 (d) part of the hyrbid
function Hsdp1=HsdtoQp1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=0 (s part of the hybrid)
Hsdp1(:,:,1)=(-1*sqrt(3)/(16*pi))*(1/(Q*r))*(4*sin((Sq/2)*r)+((pi^2-Q*r^2)/pi)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
%L=2 (d part of the hybrid)
Hsdp1(:,:,2)=(-1*sqrt(6)/(64*pi^2))*(1/(Sq^3*r^3))*(48*pi*Sq*r*cos((Sq/2)*r)+4*pi*(-24+3*pi^2-Q*r^2)*sin((Sq/2)*r)+(3*pi^4-4*pi^2*Q*r^2+Q^2*r^4)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
end

%Transition between Hsd1->Q(l=2,m=1): HQsdtoD1
%To Q(l=2) contribute 2 terms (m=1)
%Sandwitch: 1st row L=0 (s) and 2nd row L=2 (d) part of the hyrbid
function Hsdd1=HsdtoQd1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hsdd1=[];
%L=0 (s part of the hybrid)
Hsdd1(:,:,1)=(-1*sqrt(15)/(16*pi))*(1/(Sq^3*r^2))*(-8+8*cos((Sq/2)*r)+4*Sq*r*sin((Sq/2)*r)+(pi^2-Q*r^2)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)-2*sinint(pi/2)));
%L=2 (d part of the hybrid)
Hsdd1(:,:,2)=(-1*sqrt(30)/(64*pi))*(1/(Sq^5*r^4))*(576-24*pi^2+32*Q*r^2+(-576+24*pi^2+40*Q*r^2)*cos((Sq/2)*r)+(-288*Sq*r+12*pi^2*Sq*r-4*Sq^3*r^3)*sin((Sq/2)*r)+(-6*pi^4+8*pi^2*Q*r^2-2*Q^2*r^4)*sinint(pi/2)+(3*pi^4-4*pi^2*Q*r^2+Q^2*r^4)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)));
end

%Transition between Hsd1->Q(l=2,m=2): HQsdtoD2
%To Q(l=2) contribute 2 terms (m=2)
%Sandwitch: 1st row L=0 (s) and 2nd row L=2 (d) part of the hyrbid
function Hsdd2=HsdtoQd2trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=0 (s part of the hybrid)
Hsdd2(:,:,1)=0*r;
%L=2 (d part of the hybrid)
Hsdd2(:,:,2)=(1*3*sqrt(15)/(64*pi))*(1/(Sq^5*r^4))*(192-8*pi^2+16*Q*r^2+8*(-24+pi^2+Q*r^2)*cos((Sq/2)*r)-4*Sq*r*(24-pi^2+Q*r^2)*sin((Sq/2)*r)+(-2*pi^4+4*pi^2*Q*r^2-2*Q^2*r^4)*sinint(pi/2)+(pi^2-Q*r^2)^2*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)));
end


%%%%%%%%%% d2
%Transition between Hd2->Q(l=0): HQd2toS0
%Only contribution for Q(l=0)
%Sandwitch directe
function Hd2s0=Hd2toQs0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hd2s0=(1*sqrt(5)/(16*pi))*(1/(Sq^3*r^2))*(-8+8*cos((Sq/2)*r)+4*Sq*r*sin((Sq/2)*r)+(pi^2-Q*r^2)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)-2*sinint(pi/2)));
end

%Transition between Hd2->Q(l=1): HQd2toP0
%To Q(l=1) contribute 2 terms (m=0)
%Sandwitch directe
function Hd2p0=Hd2toQp0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hd2p0=(1*sqrt(15)/(16*pi))*(1/(Q^2*r^3))*(16*pi*Sq*r*cos((Sq/2)*r)+4*(-8+pi^2)*sin((Sq/2)*r)+pi*(pi^2-Q*r^2)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
end

%Transition between Hd2->Q(l=1): HQd2toP1
%To Q(l=1) contribute 2 terms (m=1)
%Sandwitch directe
function Hd2p1=Hd2toQp1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hd2p1=(-1*sqrt(6)/(64*pi^2))*(1/(Q^2*r^3))*(48*pi*Sq*r*cos((Sq/2)*r)+4*pi*(-24+3*pi^2-Q*r^2)*sin((Sq/2)*r)+(3*pi^4-4*pi^2*Q*r^2+Q^2*r^4)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
end

%Transition between Hd2->Q(l=2): HQd2toD0
%To Q(l=2) contribute 3 terms (m=0)
%Sandwitch directe
function Hd2d0=Hd2toQd0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hd2d0=(1*sqrt(5)/(32*pi))*(1/(Sq^5*r^4))*(576-24*pi^2+32*Q*r^2+(-576+24*pi^2+40*Q*r^2)*cos((Sq/2)*r)+(-288*Sq*r+12*pi^2*Sq*r-4*Sq^3*r^3)*sin((Sq/2)*r)+(-6*pi^4+8*pi^2*Q*r^2-2*Q^2*r^4)*sinint(pi/2)+(3*pi^4-4*pi^2*Q*r^2+Q^2*r^4)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)));
end

%Transition between Hd2->Q(l=2): HQd2toD1
%To Q(l=2) contribute 3 terms (m=1)
%Sandwitch directe
function Hd2d1=Hd2toQd1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hd2d1=(-1*sqrt(30)/(64*pi^2))*(1/(Q^2*r^3))*(576-24*pi^2+32*Q*r^2+(-576+24*pi^2+40*Q*r^2)*cos((Sq/2)*r)+(-288*Sq*r+12*pi^2*Sq*r-4*Sq^3*r^3)*sin((Sq/2)*r)+(-6*pi^4+8*pi^2*Q*r^2-2*Q^2*r^4)*sinint(pi/2)+(3*pi^4-4*pi^2*Q*r^2+Q^2*r^4)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)));
end

%Transition between Hd2->Q(l=2): HQd2toD2
%To Q(l=2) contribute 3 terms (m=2)
%Sandwitch directe
function Hd2d2=Hd2toQd2trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
Hd2d2=(1*3*sqrt(15)/(64*pi))*(1/(Sq^5*r^4))*(192-8*pi^2+16*Q*r^2+8*(-24+pi^2+Q*r^2)*cos((Sq/2)*r)-4*Sq*r*(24-pi^2+Q*r^2)*sin((Sq/2)*r)+(-2*pi^4+4*pi^2*Q*r^2-2*Q^2*r^4)*sinint(pi/2)+(pi^2-Q*r^2)^2*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)));
end


%%%%%%%%%% (p/f)2
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
Hpfs0(:,:,1)=(1*sqrt(2)/(32*pi))*(1/(Q*r))*(4*sin((Sq/2)*r)+((pi^2-Q*r^2)/pi)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
%L=3 (f part of the hybrid)
Hpfs0(:,:,2)=(1*sqrt(2)/(64*pi^2))*(1/(Q^2*r^3))*(80*pi*Sq*r*cos((Sq/2)*r)+4*pi*(-40+5*pi^2-Q*r^2)*sin((Sq/2)*r)+(5*pi^4-6*pi^2*Q*r^2+Q^2*r^4)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
end

%Transition between Hpf2->Q(l=1,m=0): HQpftoP0
%To Q(l=1) contribute 2 terms (m=0)
%Sandwitch: 1st row L=1 (p) and 2nd row L=3 (f) part of the hyrbid
function Hpfp0=HpftoQp0trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=1 (p part of the hybrid)
Hpfp0(:,:,1)=(1*sqrt(6)/(32*pi))*(1/(Sq^3*r^2))*(-8+8*cos((Sq/2)*r)+4*Sq*r*sin((Sq/2)*r)+(pi^2-Q*r^2)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)-2*sinint(pi/2)));
%L=3 (f part of the hybrid)
Hpfp0(:,:,2)=(1*sqrt(6)/(64*pi))*(1/(Sq^5*r^4))*(960-40*pi^2+48*Q*r^2+8*(-120+5*pi^2+9*Q*r^2)*cos((Sq/2)*r)+(-480*Sq*r+20*pi^2*Sq*r-4*Sq^3*r^3)*sin((Sq/2)*r)+(-10*pi^4+12*pi^2*Q*r^2-2*Q^2*r^4)*sinint(pi/2)+(5*pi^4-6*pi^2*Q*r^2+Q^2*r^4)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)));
end

%Transition between Hpf2->Q(l=1,m=1): HQpftoP1
%To Q(l=1) contribute 2 terms (m=1)
%Sandwitch: 1st row L=1 (p) and 2nd row L=3 (f) part of the hyrbid
function Hpfp1=HpftoQp1trans(x,Ei,Ef,M)
% x is the system = r in the string
% E=Ef-Ei of the trensition
% M = dipion invariant mass
E=Ef-Ei;
Sq=sqrt(E^2-M^2);
Q=E^2-M^2;
r=diag(x);
%L=1 (p part of the hybrid)
Hpfp1(:,:,1)=(-1*3*sqrt(2)/(32*pi))*(1/(Sq^3*r^2))*(-8+8*cos((Sq/2)*r)+4*Sq*r*sin((Sq/2)*r)+(pi^2-Q*r^2)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)-2*sinint(pi/2)));
%L=3 (f part of the hybrid)
Hpfp1(:,:,2)=(-1*sqrt(3)/(32*pi))*(1/(Sq^5*r^4))*(960-40*pi^2+64*Q*r^2+8*(-120+5*pi^2+7*Q*r^2)*cos((Sq/2)*r)+(-480*Sq*r+20*pi^2*Sq*r-12*Sq^3*r^3)*sin((Sq/2)*r)+(-10*pi^4+16*pi^2*Q*r^2-6*Q^2*r^4)*sinint(pi/2)+(5*pi^4-8*pi^2*Q*r^2+3*Q^2*r^4)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)));
end

%Transition between Hpf2->Q(l=2,m=0): HQpftoD0
%To Q(l=2) contribute 3 terms (m=0)
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
Hpfd0(:,:,1)=(1*sqrt(10)/(64*pi^2))*(1/(Q^2*r^3))*(48*pi*Sq*r*cos((Sq/2)*r)+4*pi*(-24+3*pi^2-Q*r^2)*sin((Sq/2)*r)+(3*pi^4-4*pi^2*Q*r^2+Q^2*r^4)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
%L=3 (f part of the hybrid)
Hpfd0(:,:,2)=(1*sqrt(10)/(128*pi^2))*(1/(Q^3*r^5))*(16*pi*Sq*r*(-720+15*pi^2+7*Q*r^2)*cos((Sq/2)*r)+4*pi*(15*(384-8*pi^2+pi^4)-8*(67+pi^2)*Q*r^2+Q^2*r^4)*sin((Sq/2)*r)+(15*pi^6-23*pi^4*Q*r^2+9*pi^2*Q^2*r^4-Q^3*r^6)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
end

%Transition between Hpf2->Q(l=2,m=1): HQpftoD1
%To Q(l=2) contribute 3 terms (m=1)
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
Hpfd1(:,:,1)=(-1*3*sqrt(10)/(64*pi^2))*(1/(Q^2*r^3))*(16*pi*Sq*r*cos((Sq/2)*r)+4*(-8+pi^2)*sin((Sq/2)*r)+pi*(pi^2-Q*r^2)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
%L=3 (f part of the hybrid)
Hpfd1(:,:,2)=(-1*sqrt(15)/(32*pi))*(1/(Q^3*r^5))*(16*Sq*r*(-240+50*pi^2+2*Q*r^2)*cos((Sq/2)*r)+(20*(384-8*pi^2+pi^4)-4*(176+3*pi^2)*Q*r^2)*sin((Sq/2)*r)+pi*(5*pi^4-8*pi^2*Q*r^2+3*Q^2*r^4)*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
end

%Transition between Hpf2->Q(l=2,m=2): HQpftoD2
%To Q(l=2) contribute 3 terms (m=2)
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
Hpfd2(:,:,1)=(1*sqrt(15)/(64*pi^2))*(1/(Q^2*r^3))*(16*pi*Sq*r*cos((Sq/2)*r)+4*pi*(-8+pi^2-Q*r^2)*sin((Sq/2)*r)+(pi^2-Q*r^2)^2*(sinint((pi-Sq*r)/2)-sinint((pi+Sq*r)/2)));
%L=3 (f part of the hybrid)
Hpfd2(:,:,2)=(-1*sqrt(10)/(128*pi^2))*(1/(Q^3*r^5))*(16*pi*Sq*r*(240-5*pi^2+Q*r^2)*cos((Sq/2)*r)+-4*pi*(1920+5*pi^4-152*Q*r^2+Q^2*r^4-2*pi*(20+3*Q*r^2))*sin((Sq/2)*r)-(pi^2-Q*r^2)^2*(5*pi^2-Q*r^2)*(sinint((pi-Sq*r)/2)+sinint((pi+Sq*r)/2)));
end
