% The function GGraph is generating a directed or undirecte Bernoulli graph 
% which has community structures via Stochastic Blockmodels.
% Input: n is the number of nodes.
%        Pi is the propotion of the different groups.
%        P is the probability matrix.
%        gtype is "directed" or "undirected".

% Output: A is the n-by-n adjacent matrix of a graph.
%         Glabel is the group label of nodes.
%         Z is the indicator of the group of nodes.
function [A, Glabel, Z] = GGraph(n, Pi, P, gtype)

A = zeros(n);

Glabel = zeros(1,n);

Z = mnrnd(1, Pi, n);% get indicators

[temp1, temp2] = find(Z);

Glabel(temp1) = temp2; % get group labels

switch gtype
    case 'directed'
        A = (Z*P*Z'<=rand(n));
        A = A - diag(diag(A));
    case 'undirected'
        A = (Z*P*Z'<=rand(n));
        A = triu(A) - diag(diag(A));
        A = A + A';
end




