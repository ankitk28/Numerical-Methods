 function [w, x] = gaussint(n)
 format long
     if (n==1)
         x = 0;
         w = 2;
     else
         xcurr = 0;
         x = zeros(n,1);
         for i = 2:n
             xnext = zeros(i,1);
             for j = 1:i-1
                 if(j == 1)
                     a = -1;
                     b = xcurr(1);
                 else
                     a = xcurr(j-1);
                     b = xcurr(j);
                 end
                 xnext(j) = schroderbisection(a, b, i, 0.000000000001);
             end
             xnext(i) = schroderbisection(xcurr(i-1), 1, i, 0.000000000001);
             xcurr = xnext;
             if(i==n)
                 x = xcurr;
             end
         end
         w = zeros(n,1);
         for i = 1:n
            basis = @(y)lagrangesquare(y, x, i);
            w(i) = gadap(-1, 1, basis, 0, 0.000000000001);
         end
     end
     
     function val = lagrangesquare(eva_at,inp_set,order)
        len = numel(inp_set);
        li = 1;
        for j1 = 1:len
            if(order ~= j1)
                li = (li*(eva_at-inp_set(j1)))/(inp_set(order)-(inp_set(j1)));
            end
        end
        val = li*li;
     end

 end
 