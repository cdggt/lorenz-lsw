function m = orbit_mean(x,dim)
%ORBIT_MEAN This method computes the average of PERIODIC data x, over
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
m = fft(x,[],dim);
S.subs = repmat({':'},1,ndims(m));
S.subs{dim} = 1;
S.type = '()';
m = subsref(m,S)/N;

end