function [Elow,Eup,status] = getBracket(system,md,k,E,step)
% [Elow,Eup] = getBracket(system,md,k,E,step)
%  computes a bracket [elow,eup] containing the eigenvalue with index k.
%  The bracket is sufficiently small such that values in this interval form
%  good starting values for a shooting eigenvalue search


status=1;
[m,st] =getIndex(system,md,E);
status= min(status,st);
while true
  if m <= k 
      mlow = m;
	  Elow = E;
	  E = E + step;
      step = step*2;
  else
	  Eup = E;
      Elow=0;
      mup=m;
	  break;
  end
  [m,st] = getIndex(system,md,E);  status= min(status,st);
end
%now we have found Elow and Eup such that Elow <E_k< Eup


its=0;
while true
    E = (Elow + Eup)/2;
    if mup>k+3
       E = Elow+(Eup - Elow)/4;
    end
    m = getIndex(system,md,E);
    mt = getIndex(system,md,E+(md.tol)^(1/3));
    while mt ~= m
        E=E+(md.tol)^(1/3);
        m=mt;
        mt=getIndex(system,md,E+(md.tol)^(1/3));
    end
    if m <=k
      Elow=E;
      if m==k && mup==k+1 
        return;
      end
      if mlow==m 
          its=its+1;
          if its>12
             return;
          end
      end
      mlow=m;
    else
        Eup=E;
        if mup==m
            its=its+1;
            if its>10 %probably means that the multiplicity > 1
                return;
            end
        end
        mup=m;
    end
    if abs(Eup-Elow)<max((md.tol)^(1/3),1e-10)  %1/3 since h^6->tol, h^2->tol^(1/3)
        return;
    end
end

end

