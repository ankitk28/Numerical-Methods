function c = cleg(n)
format long
% n: number of recurrence parameters
if(n==0)
    c = [];
elseif( n== 1)
    c = 0;
elseif( n== 2)
    c = [0 1/3];
else
    
     c = cleg(n-1);
%      numel(c)
%      disp('k');

     pnminusnume = @(x) x.*value(x,c(1:n-1)).*value(x,c(1:n-2));
    
     numer = gadap(-1,1,pnminusnume,0,0.000000000001);
      
     pnminus2 = @(x) valuesquare(x,c(1:n-2));
     deno = gadap(-1,1,pnminus2,0,0.000000000001 );
     
     c(n) = numer/deno;  
     %disp('k');
end

function val = value(x,car)
  [val, zz, zzz] = pleg(x,car);
end

function val = valuesquare(x,car)
  [valsq, zz, zzz] = pleg(x,car);
  val = valsq*valsq;
end




end
    