function [ attn ] = getAttn(nSteps)
%GETATTN Summary of this function goes here
%   Detailed explanation goes here


% Initial attention
attn(1) = rand();

% Delay to max attention
maxdel = exprnd(5);
attn(maxdel) = 1;

% Min attention



end

