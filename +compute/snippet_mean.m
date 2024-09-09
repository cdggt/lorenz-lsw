function m = snippet_mean(x,dim)
%SNIPPET_MEAN This method computes the average of NONPERIODIC data x, over
%dimension dim
%
% Inputs:
%
%   x    : an N-D array
%   dim  : (optional, default = 1) an integer. the dimension along which to 
%        compute mean

if nargin < 2
    dim = 1;
end

N = size(x,dim);
m = trapz(x,dim)/(N-1);

end