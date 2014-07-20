% Generate a graph
function [A, Glabel,Z] = GGraph(n,Alpha,PI,Distribution,NetType)

A = zeros(n);
Z = mnrnd(1,Alpha,n);

Glabel = zeros(1,n);
for i = 1:n
    Glabel(i) = find(Z(i,:)==1);
end

switch Distribution
    
    case 'Bernoulli'
        
        switch NetType
            
            case 'directed'
                for i = 1:n
                    for j = 1:n
                        edge_seed = rand(1);
                        if edge_seed<= PI(Glabel(i),Glabel(j))
                            A(i,j) = 1;
                        end
                    end
                end
                
            case 'undirected'
                for i = 2:n
                    for j = 1:(i-1)
                        edge_seed = rand(1);
                        if edge_seed<= PI(Glabel(i),Glabel(j))
                            A(i,j) = 1;
                        end
                    end
                end
        end
        A = A - diag(diag(A));
        A = A + A';
        
    case 'Poisson'
        
        switch NetType
            
            case 'directed'
                for i = 1:n
                    for j = 1:n
                        x = 1:15;
                        pro = poisspdf(x,PI(Glabel(i),Glabel(j)));
                        [value,index] = sort(pro,'descend');
                        edge_seed = rand(1);
                         EdgeNumber = cal_edge_num(value,index,edge_seed);
                        A(i,j) = EdgeNumber;
                    end
                end
                A  = A - diag(diag(A));
            case 'undirected'
                for i = 2:n
                    for j = 1:i-1
                        x = 1:15;
                        pro = poisscdf(x,PI(Glabel(i),Glabel(j)));
                        [value,index] = sort(pro,'ascend');
                        tmp = exp(-PI(Glabel(i),Glabel(j)));
                        edge_seed = rand(1);
                         EdgeNumber = cal_edge_num(tmp,value,index,edge_seed);
                        A(i,j) = EdgeNumber;
                    end
                end
                A  = A - diag(diag(A));
                A = A + A';
        end
end

% B = Z * PI * Z';
% C = rand(n);
% A = (C<=B);
% A = A - diag(diag(A));