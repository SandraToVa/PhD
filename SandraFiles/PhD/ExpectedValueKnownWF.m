
function V=ExpectedValueKnownWF(Wi,Wf,x,Op,s)
%This code is the same for Expected value but in this case we know the wave
%functions as a vector of x dimension (must be the same for the 2 wf)

%Function to compute espected values using wavefunctions and any observable
%Output:        - Returns a vector of positions = columns of wavefunctions initial
%               and final (if the hybrid computed is a 3x3 matrix, 3 colums, if it is a
%               7x7 matrix, 7 columns)
%               This vector is the expectation value of
%               final wavefunction * operator * initial wavefunction
%Input:         - Wi = wave function of the initial state (a matrix with
%               the proper length "xi"
%               - Wf= wave function of the final state (a matrix with
%               the proper length "xf"=xi
%               - x= the lenght of the wf =xi=xf
%               - Op= operator in matrix form (it has to be a diagonal
%               matrix) the rowas and columns the same as the
%               wavefunctions. This is the operator that we want to find
%               the expectation value of
%               - s= spin

%I have to modify the tolerance in "initial" and "final" for the mesh to
%coincide. And for the x and lenght of the 2nd dimension of W to coincide!

%The length of the wavefunctions must be = x
[numRi,~] = size(Wi);
[numRf,~] = size(Wf);

if numRi ~= numRf
    disp('Number of rows in wave functions not equal')
end

%If we work with quarkonium the wave functions only have on relevant row:
%the first one
if s==0
    Wi_f=Wi(1,:); %row vector
    Wf_f=Wf(1,:); %row vector
    
    V1=Wf_f*Op;
    V2=zeros(1,length(V1)); %row of zeros
    for el=1:length(V1)
        %Second auxiliar vector a row vector. Final multiplication of the
        %Yf*Op*Yi
        V2(el)=V1(el)*Wi_f(el);
    end
    %The integrand is exactly V2
    I=V2;
    V=trapz(x,I);
     
end %end of if
        
%If I work with hybrid-to-qurakonium, the hybrid wave functions have 2 rows
%for each row I have to apply different operators between the quarkonium
%and the bybrid wavefunctions

if s==1
    %If the initial states is hybrid, (should be) numRi=3 if the hybrid is
    %a double the two first rows are treated separately and if the hybrid
    %is not a double the only relevant row is the last one
    %If the final state is quarkonium, (should be) numRi=3 but we only
    %want the first one
        
    V=zeros(numRi,1);
      
    for N=1:numRi
        %Muliply each row for any row of the other wavefunction normalized
        
        Wi_f=Wi(N,:); %row vector
        Wf_f=Wf(1,:); %row vector

        %The operator Op in this case is three elements, one for the upper row and
        %one for the bottom row (in the case of double) and the third for
        %the case of single
        OpRow=Op(:,:,N);

        %With out wavefunctions computed now we multiply for our operator
        %(matrix) Yf*Op*Yi
        %Auxiliar matrix(vector) for the first multiplication. It is a row
        %vector
        V1=Wf_f*OpRow;
        %Now our second multiplication is for each element in V1 and each
        %element in Yi
        V2=zeros(1,length(V1)); %row of zeros
        for el=1:length(V1)
           %Second auxiliar vector a row vector. Final multiplication of the
           %Yf*Op*Yi
           V2(el)=V1(el)*Wi_f(el);
        end
        %The integrand is exactly V2
        I=V2;
        %We obtain a column resulting of the sandwitch of
        %wavefunctionQ*Op[N]*wavefunctionH[N]

        V(N,1)=trapz(x,I);  
        
    end %end for

    %Then we should add the two results. But we do not do this here
    
end %end if

end %end function


