 function [Est_Glabel,Tau,Alpha,W,Theta,iter,time] = VEM_ODSBM(A,Q,maxIter,tol,init_type)
% Usage: Variational EM algorithm for DegreeCorrected Overlapping Stochastic Block
% Model
% input
% A - adjacent matrix
% Q -  the number of grounps
% maxIter - maximum number of iterations
%         - DEFAULT 500
% tol - tolerance for stopping criterion
%     - DEFAULT 1e-6
% output
% Est_Glabel: the group labels of each vertice
% Tau - n x Q
%     - Tau_iq  -  the probability of vertice i belongs to class q
% Alpha - 1 x Q
%    - Alpha_q  -  the fraction of vertices in class q
% W  - Q x Q
%    - W_ql  - the average number of edges between a vertice in class q 
%            - and a vertice in class l
% Theta - 1 x n
%       - Theta_i -  parameter to control the degree of the i_th vertice
% iter - the number of iteration
% time - time cost of the algorithm

if nargin < 5
    init_type = 'spectral';
end
if nargin < 4
    tol = 1e-6;
end
if nargin < 3
    maxIter = 500;
end

n = size(A,1);   % the number of vertices

%initialization Tau, Theta, W
switch init_type
    case 'random'
        tmp = rand(n,Q);
       Tau = tmp./repmat(sum(tmp,2),1,Q);
    case 'spectral'
        [ C, ~, ~ ] = SpectralClustering( A, Q, 3 );
        Tau = full(C);
end
temp = rand(1,n);
Theta = temp/sum(temp);
W = rand(Q);

iter = 1;
while  iter<=maxIter
    iter
    %----------- M step -------------
    %update Alpha
    Alpha = sum(Tau)/n;
    %update Theta 
    in_iter = 1;
    int_tol = 1e-4;
    diff = int_tol + 1;
    max_InIter = 10;
    while in_iter< max_InIter && diff>int_tol
        old_Theta = Theta;
        tmp1 = sum(A,2);
        tp1 = Tau*W*Tau';
        tp2 = repmat(Theta,n,1);
        tmp2 = sum(tp2.*tp1,2);
        Theta = tmp1./tmp2;
        diff = max(max(abs(Theta-old_Theta)));
        in_iter = in_iter + 1;
    end
    
    %update W
    W = Gradient_Ascent (W,A,Tau,Theta);
    
     % convergence
    if iter > 1
        old_L = L;
    end
    L = LBLogComlikelihood(A,Alpha,Theta,W,Tau);
    if iter > 1
        err = abs(1-L/old_L);
        if err < tol
            break;
        end
    end
    
    %update Tau
    
    
    
    iter = iter + 1;
end
[~,Est_Glabel] = max(Gamma,[],2);
