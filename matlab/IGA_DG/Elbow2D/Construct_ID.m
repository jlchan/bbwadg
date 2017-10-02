%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function: Construct_ID
%
% Input:  d = number of spatial dimensions
%         n = number of functions
%
% Output: ID = ID array
%
% Purpose: Construct the ID array
%
% Notes: N/A
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ID] = Construct_ID(d,n)

for A = 1:d    
    for i = 1:n
        ID(A,i) = (i-1)*d+A;
    end
end