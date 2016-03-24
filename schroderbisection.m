function [r, h] = schroderbisection(a, b, n, t)
format long
% Function that combines the fast convergence of Schroder iteration for
% multiple roots with the bracketing guarantee of bisection
%
% a: beginning of interval [a, b]
% b: end of interval [a, b]
% f: function handle f(x) to find a zero for
% fp: function hanfle f'(x)
% fpp: function handle f''(x)
% t: user provided tolerance for interval width

% the left endpoint of our interval
left = a;
% the right endpoint of our interval
right = b;
% the size of the interval
range = right - left;

% the matrix within which we will store the various values we acquire
% as we progress through the iterations
h = zeros(6, 1);

% a counter to keep track of the number of iterations
count = 1;
% the function e as defined in the specification
car = cleg(n);
e = @(x) min( abs( value(right,car)  - value(left,car) ) /8 , abs( dder(x,car) ) * abs(right - left)^2);

% we create two separate iteration functions - one for g+ and another for
% g- which are simply the regular iteration function except with f(x)
% replaced with (f(x) + e(p)) in the g+ case and (f(x) - e(p)) in the
% g- case.
gplus = @(x) x - (((value(x,car) + e(x)) * der(x,car)) / (der(x,car)^2 - ((value(x,car) + e(x)) * dder(x,car))));
gminus = @(x) x - (((value(x,car) - e(x)) * der(x,car)) / (der(x,car)^2 - ((value(x,car) - e(x)) * dder(x,car))));

while right - left >= t && right - left >= eps(left)    
    h(1, count) = left;
    h(5, count) = right;
    
    % we define p at each step of the iteration through geometric
    % bisection
    if left <= 0 && right >= 0
        if sign(value(realmin,car)) ~= sign(value(right,car))
            p = realmin;
        else
            p = -realmin;
        end
    else
        p = sqrt(abs(left * right)) * sign(left);
    end
    
    h(2, count) = p;
    h(3, count) = gminus(p);
    h(4, count) = gplus(p);
    h(6, count) = value(p,car);
    
    % the five values from which we will select our new interval. We then
    % sort it and add it into a new list.
    values = [left, p, gminus(p), gplus(p), right];
    sorted = sort(values);
    
    % if f(v) = 0 for any of the five values above, then we have found
    % the root and set both the left and right endpoint to be this value
    for i = 1:5
        if value(sorted(i),car) == 0
            left = sorted(i);
            right = sorted(i);
            break
        end
    end
    
    % if the left and right endpoints are identical, then we have found
    % the root and we can exit the loop
    if left == right
        break
    end
    
    % we now iterate through adjacent pairs of values and select the first 
    % pair where we observe a sign change - this interval will therefore
    % contain the root and since this list of values is sorted, we know
    % that this is the smallest interval containing the interval. We choose
    % to sort the list and iterate through it this way instead of iterating
    % through every pair and picking the smallest interval because if any
    % of the values are realmin or smaller than realmin, MATLAB will not be
    % able to distinguish between the two values and we could end up with a
    % bug where the code cannot distinguish between the two intervals.
    for i = 1:4
        if sign(value(sorted(i),car)) ~= sign(value(sorted(i + 1),car))
            left = sorted(i);
            right = sorted(i + 1);
            break;
        end
    end
    
    count = count + 1;
end

% now that we have found the required interval, we store the endpoints of
% the interval and the other required values in our matrix
if left <= 0 && right >= 0
    p = 0;
else
    p = sqrt(abs(left * right)) * sign(left);
end

h(1, count) = left;
h(2, count) = p;
h(3, count) = gminus(p);
h(4, count) = gplus(p);
h(5, count) = right;
h(6, count) = value(p,car);

% we compute the approximate to the root, r, as the arithmetic mean of
% our left and right endpoint
r = (left + right)/2;


function val = value(x,car)
  [val, zz, zzz] = pleg(x,car);
end

function vald = der(x,car)
  [z, vald, zzz] = pleg(x,car);
end

function valdd = dder(x,car)
  [z, zz, valdd] = pleg(x,car);
end


end
