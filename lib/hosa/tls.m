function [x,flag] = tls(A,b)
%TLS	Total-Least Squares solution to a set of linear equations
%	[x,flag] = tls (A, b)
%	Obtains the TLS solution to  A x = b
%	Must have: # rows (A) >= # cols (A) + # cols(b)
%	flag is set to 1, if the problem does not have a unique solution.

%  Copyright (c) 1991-2001 by United Signals & Systems, Inc. 
%       $Revision: 1.4 $
%  A. Swami   January 20, 1993.

%     RESTRICTED RIGHTS LEGEND
% Use, duplication, or disclosure by the Government is subject to
% restrictions as set forth in subparagraph (c) (1) (ii) of the
% Rights in Technical Data and Computer Software clause of DFARS
% 252.227-7013.
% Manufacturer: United Signals & Systems, Inc., P.O. Box 2374,
% Culver City, California 90231.
%
%  This material may be reproduced by or for the U.S. Government pursuant
%  to the copyright license under the clause at DFARS 252.227-7013.

[mrow,ncol]=size(A);
[nrow,kcol]=size(b);

if (mrow ~= nrow)
    error('A, b should have same number of rows')
end

if (mrow < ncol + kcol)
   error('must have row-dim(A) >= col-dim([A,b]) ')
end

[u,s,v] = svd([A,b]);  s = diag(s);
if (s(ncol) == s(ncol+1))
    disp('TLS solution may be singular')
    flag = 1;
end

v12 = v(1:ncol,ncol+1:ncol+kcol);
v22 = v(ncol+1:ncol+kcol,ncol+1:ncol+kcol);
x   = -v12/v22;

return
