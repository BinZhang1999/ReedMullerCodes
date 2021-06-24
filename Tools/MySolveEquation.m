function [A,x] = MySolveEquation(A, b)
%% Solve equation Ax=b over gf(2)
% Input: 
%   A: mXn matrix. rank(A)=n
%   b: nX1 vector.
% Output:
%   x: mX1 vector
%%
[m,n]=size(A);

for i = 1:n
   idx = find(A(i:end,i)==1);
   if isempty(idx)
       break;
   end
   idx = idx+i-1;
   
   vec_x = A(idx(1),:); b_x = b(idx(1));
   A(idx(1),:) = A(i,:); b(idx(1)) = b(i);
   A(i,:) = vec_x; b(i) = b_x;
   
   for j = 1:m
       if j==i
           continue;
       end
       if (A(j,i)==1)
           A(j,:)=(A(j,:)~=A(i,:)); b(j) = (b(j)~=b(i));
       end
   end
end
x = b;
end