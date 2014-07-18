function [PI, Alpha, Tau, CluResult,time] = VEM(X, Q, Distribution,MaxIter,IniType, NetType, )
%
% Author - Xiangyong Cao, 07/2014
%
% Email  - caoxiangyong45@gmail.com
%
% Description - main function to implement the variational EM algorithm.
%
% Input  - X    : adjacency matrix   n x n
%        - Q    : cluster number
%        - MaxIter   : iteration number, default 200
%        - IniType   : 'random','spectral'
%        - NetType   : 'undirected','directed'
%        - Distribution: 'Bernoulli','Possion'
% Output - PI        : for 'Bernoulli', the link probability matrix   Q x Q
%                    : for 'Possion', the average number of edges between
%                    groups
%        - Alpha     : proportion of each group
%        - Tau       : posterior probability of each vertex
%        - CluResult : label of each vertex
%        - time       : optimizing time
% -------------------------------------------------------------------------

% default settings-----------------------------


if nargin < 6    
    NetType = 'undirected';
end
if nargin < 5    
    IniType = 'Spectral';
end
if nargin < 4    
    MaxIter = 200;           
end
if nargin < 3    
    Distribution = 'Bernoulli';           
end

% common constants---------------------
n = size(X,1);   

% Init --------------------------------------------
switch IniType
    case 'random'
        tmp = rand(n,Q);
        Tau = tmp./repmat(sum(tmp,2),1,Q);
    case 'spectral'
        [ C, L, U ] = SpectralClustering( X, Q, 3 );
        Tau = full(C);
end

% ----------internal control parameters----------
int_maxIter = 5;
out_tol = 1e-10;
int_tol = 1e-4;
log_min = log(realmin);
prune   = 1e-10;


tic;
out_iter = 1;
while out_iter<MaxIter
    
    out_iter
    
    %% Variational M step
    % update Alpha, PI, Theta      
    col_sum = sum(Tau);
    Alpha = col_sum/n;
    PI = (Tau' * X * Tau)./(col_sum'*col_sum - Tau'*Tau);  
    
    % convergence
        if out_iter > 1       
            old_L = L; 
        end
        L = lblikelihood(X,Tau,Alpha,PI,NetType);
        if out_iter > 1
            err = abs(1-L/old_L);                                                               % relative error of lower bound
            if err < out_tol 
                break;
            end
        end            
        
    %% Variational E step
    % updata Tau
    A = repmat(log(Alpha),n,1);
    b = log(1-PI);
    B1 = log(PI) - b;
    B1(B1==-inf) = 0;
    int_iter = 1;
    diff = int_tol + 1;
    while (int_iter<int_maxIter) && (diff>int_tol) &&(keep_on)
        old_Tau = Tau;
        tmp = Tau*b';
        tmp = repmat(sum(tmp,1),n,1)-tmp;
        tmp1 = X*Tau*B1';
        switch NetType
            case 'undirected'
                Tau = 0.5*(tmp1+tmp)+A;
            case 'directed'
                Tau= tmp1+tmp+A;
        end
        Tau(Tau<log_min) = log_min;
        Tau = Tau-repmat(max(Tau,[],2),1,Q);
        Tau = exp(Tau);
        Tau = Tau./repmat(sum(Tau,2),1,Q);
        Tau(Tau<prune) = 0;
        if int_iter > 1
            old_diff = diff;
            diff = max(max(abs(Tau-old_Tau)));
            if diff >= old_diff
                keep_on = false;                                                                                % diff increase, stop iterating
            end
        elseif int_iter == 1
            diff = max(max(abs(Tau-old_Tau)));
        end
        int_iter = int_iter + 1;
    end
        out_iter = out_iter + 1;
end
time = toc;
%CluResult
if Q == 1
    CluResult = ones(1,n);
else
    CluResult = CluResult1;
end


