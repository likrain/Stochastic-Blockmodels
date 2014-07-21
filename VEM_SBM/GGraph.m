% The function GGraph is generating a directed or undirecte Bernoulli and Poisson graph
% which has community structures via Stochastic Blockmodels.
% Input: n is the number of nodes.
%        Alpha is the propotion of the different groups.
%        PI is the probability matrix.
%        Distribution is "Bernoulli" or "Poisson"
%        NetType is "directed" or "undirected".

% Output: A is the n-by-n adjacent matrix of a graph.
%         Glabel is the group label of nodes.
%         Z is the indicator of the group of nodes.
function [A, Glabel,Z] = GGraph(n,Alpha,PI,Distribution,NetType)

A = zeros(n);
Z = mnrnd(1,Alpha,n);

Glabel = zeros(1,n);
[temp1, temp2] = find(Z);
Glabel(temp1) = temp2; % get group labels

switch Distribution
    case 'Bernoulli'
        switch NetType
            case 'directed'
                A = (Z*PI*Z'>=rand(n));
                A = A - diag(diag(A));
            case 'undirected'
                A = (Z*PI*Z'>=rand(n));
                A = triu(A) - diag(diag(A));
                A = A + A';
        end
    case 'Poisson'
        switch NetType
            case 'directed'
                for i = 1:n
                    for j = 1:n
                        x = 1:50;
                        pro = poisscdf(x,PI(Glabel(i),Glabel(j)));
                        tmp = exp(-PI(Glabel(i),Glabel(j)));
                        edge_seed = rand(1);
                        EdgeNumber = cal_edge_num(tmp,pro,x,edge_seed);
                        A(i,j) = EdgeNumber;
                    end
                end
                A  = A - diag(diag(A));
            case 'undirected'
                for i = 2:n
                    for j = 1:i-1
                        x = 1:50;
                        pro = poisscdf(x,PI(Glabel(i),Glabel(j)));
                        tmp = exp(-PI(Glabel(i),Glabel(j)));
                        edge_seed = rand(1);
                        EdgeNumber = cal_edge_num(tmp,pro,x,edge_seed);
                        A(i,j) = EdgeNumber;
                    end
                end
                A  = A - diag(diag(A));
                A = A + A';
        end
end