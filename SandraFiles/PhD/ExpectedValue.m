
function V=ExpectedValue(i,f,Op,initial,final,m,sin,sfin)
%Function to compute espected values using wavefunctions and any observable
%Output:        - Returns a vector of positions = columns of wavefunctions initial
%               and final (if the hybrid computed is a 3x3 matrix, 3 colums, if it is a
%               7x7 matrix, 7 columns)
%               This vector is the expectation value of
%               final wavefunction * operator * initial wavefunction
%Input:         - i=position in the spectum of eneries of the hybrid or
%               quarkonium that we want the inital wavefunction to have
%               - f= same as i but for final wavefucntion
%               - Op= operator in matrix form (it has to be a diagonal
%               matrix) the rowas and columns the same as the
%               wavefunctions. This is the operator that we want to find
%               the expectation value of
%               - initial= initial wavefunction form the scripts of "Spin..."
%               for hybrid and "Quarkonium..." for quarkonium
%               - final= same as initial for final state
%               - m= mass of the quark. Charm or Bottom

%I have to modify the tolerance in "initial" and "final" for the mesh to
%coincide. And for the x and lenght of the 2nd dimension of W to coincide!
setspin(sin)
[Ei,Wi,~]=initial(m,spin);

setspin(sfin)
[Ef,Wf,x]=final(m,spin);

%The vale of our initial state
%disp('Initial state i')
%i;
%disp('With energy:')
Ei(i);
%disp('Wavefunction:')
Yi=Wi(:,:,i);

%The length of the wavefunction i 
[numRi,~] = size(Yi);


%The vale of our final state
%disp('Initial final f')
%f;
%disp('With energy:')
Ef(f);
%disp('Wavefunction:')
Yf=Wf(:,:,f);

%The length of the wavefunction f 
[numRf,~] = size(Yf);

%If we work with quarkonium the wave functions only have on relevant row:
%the first one
if sin==0 && sfin ==0
    Yi_f=Yi(1,:); %row vector
    Yf_f=Yf(1,:); %row vector
    V1=Yf_f*Op;
    V2=zeros(1,length(V1)); %row of zeros
    for el=1:length(V1)
        %Second auxiliar vector a row vector. Final multiplication of the
        %Yf*Op*Yi
        V2(el)=V1(el)*Yi_f(el);
    end
    %The integrand is exactly V2
    I=V2;
    V=trapz(x,I);
     
end %end of if
        
%If I work with hybrid-to-qurakonium, the hybrid wave functions have 2 rows
%for each row I have to apply different operators between the quarkonium
%and the bybrid wavefunctions

if sin==1 && sfin==0
    %If the initial states is hybrid, (should be) numRi=3 if the hybrid is
    %a double the two first rows are treated separately and if the hybrid
    %is not a double the only relevant row is the last one
    %If the final state is quarkonium, (should be) numRi=3 but we only
    %want the first one
        
    V=zeros(numRi,1);
      
    for N=1:numRi
        %Muliply each row for any row of the other wavefunction normalized
        
        Yi_f=Yi(N,:); %row vector
        Yf_f=Yf(1,:); %row vector

        %The operator Op in this case is three elements, one for the upper row and
        %one for the bottom row (in the case of double) and the third for
        %the case of single
        OpRow=Op(:,:,N);

        %With out wavefunctions computed now we multiply for our operator
        %(matrix) Yf*Op*Yi
        %Auxiliar matrix(vector) for the first multiplication. It is a row
        %vector
        V1=Yf_f*OpRow;
        %Now our second multiplication is for each element in V1 and each
        %element in Yi
        V2=zeros(1,length(V1)); %row of zeros
        for el=1:length(V1)
           %Second auxiliar vector a row vector. Final multiplication of the
           %Yf*Op*Yi
           V2(el)=V1(el)*Yi_f(el);
        end
        %The integrand is exactly V2
        I=V2;
        %We obtain a column resulting of the sandwitch of
        %wavefunctionQ*Op[N]*wavefunctionH[N]

        V(N,1)=trapz(x,I);  
        
    end %end for

    %Then we should add the two results. But we do not do this here
    
end %end if

%SHA DE MIRSAR EN
%CONTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%The same if initial is quarkonium and final is hyrbid

if sin==0 && sfin==1
    %If the initial states is quarkonium, (should be) numRi=3 but only the
    %first row important
    %If the final state is hybrid, (should be) numRi=2
        
    V=zeros(numRf,1);
      
    for M=1:numRf
        %Muliply each row for any row of the other wavefunction normalized
        
        Yi_f=Yi(1,:); %row vector
        Yf_f=Yf(M,:); %row vector

        %The operator Op in this case is two elements, one for the upper row and
        %one for the bottom row.
        OpRow=Op(:,:,M);

        %With out wavefunctions computed now we multiply for our operator
        %(matrix) Yf*Op*Yi
        %Auxiliar matrix(vector) for the first multiplication. It is a row
        %vector
        V1=Yf_f*OpRow;
        %Now our second multiplication is for each element in V1 and each
        %element in Yi
        V2=zeros(1,length(V1)); %row of zeros
        for el=1:length(V1)
           %Second auxiliar vector a row vector. Final multiplication of the
           %Yf*Op*Yi
           V2(el)=V1(el)*Yi_f(el);
        end
        %The integrand is exactly V2
        I=V2;
        %We obtain a column resulting of the sandwitch of
        %wavefunctionQ*Op[N]*wavefunctionH[N]

        V(M,1)=trapz(x,I);  
        
    end %end for

    %Then we should add the two results. But we do not do this here


%Case Hybrid-Hybrid???????
    
end %end if


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

