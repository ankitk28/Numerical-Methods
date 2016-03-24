function [int,abt] = gadap(a, b, f, p, tol)
format long
abt = [];
abt(1,1)=a;
abt(1,2)=b;
tj = (   (8/9)*( f((b+a)/2)) + (5/9)*f(( b+a+sqrt(0.6)*(b-a)  )/2)   + (5/9)*f(( b+a-sqrt(0.6)*(b-a)  )/2) )*((b-a)/2);
abt(1,3) = tj;
mid = (a+b)/2;
tl = (   (8/9)*( f((mid+a)/2)) + (5/9)*f(( mid+a+sqrt(0.6)*(mid-a)  )/2)   + (5/9)*f(( mid+a-sqrt(0.6)*(mid-a)  )/2) )*((mid-a)/2); 
tr = (   (8/9)*( f((b+mid)/2)) + (5/9)*f(( b+mid+sqrt(0.6)*(b-mid)  )/2)   + (5/9)*f(( b+mid-sqrt(0.6)*(b-mid)  )/2) )*((b-mid)/2);  

if(abs(tj - (tl+tr)) > tol*max(abs(tj),abs(tl+tr)))
    [tl abtl] = gadap(a, mid, f, p, tol);
    [tr abtr] = gadap(mid, b, f, p, tol);
    ksize = size(abtl);
    for j = 1:ksize(1)
        abtsize = size(abt);
        abt(abtsize(1),:) = abtl(j,:);
    end
    ksize = size(abt);
    for j = 1:ksize(1)
        abtsize = size(abt);
        abt(abtsize(1),:) = abtr(j,:);
    end
    int = tl + tr;
else
    int = tj;
end

end
