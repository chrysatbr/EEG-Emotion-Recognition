% cite: https://www.mathworks.com/matlabcentral/answers/63247-get-max-value-and-index-of-multidimensional-array
function [row,column,value] = maxMatrix(A)
[y,~] = max(A);
[value,column] = max(y);
[~,row] = max(A(:,column));
end