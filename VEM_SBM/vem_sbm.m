function [PI, Alpha, Tau, Est_Glabel,time] = vem_sbm(X, Q, Distribution,MaxIter,IniType, NetType)
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
%        - Distribution: 'Bernoulli','Poisson'
% Output - PI        : for 'Bernoulli', the link probability matrix   Q x Q
%                    : for 'Possion', the average number of edges between
%                    groups
%        - Alpha     : proportion of each group
%        - Tau       : posterior probability of each vertex
%        - Est_Glabel : label of each vertex
%        - time       : optimizing time
% -------------------------------------------------------------------------

% default settings-----------------------------
if nargin < 6
    NetType = 'undirected';
end
if nargin < 5
    IniType = 'spectral';
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
int_maxIter = 10;
out_tol = 1e-10;
int_tol = 1e-6;
log_min = log(realmin);
prune   = 1e-10;

tic;
out_iter = 1;
while out_iter<MaxIter
    
    out_iter
    
    %% Variational M step
    % update Alpha, PI
    col_sum = sum(Tau);
    Alpha = col_sum/n;
    PI = (Tau' * X * Tau)./(col_sum'*col_sum - Tau'*Tau);
    
    % convergence
    if out_iter > 1
        old_L = L;
    end
    L = lblikelihood(X,Tau,Alpha,PI,Distribution,NetType);
    if out_iter > 1
        err = abs(1-L/old_L);
        if err < out_tol
            break;
        end
    end
    
    %% Variational E step
    % updata Tau
    switch Distribution
        case 'Bernoulli'
            A = repmat(log(Alpha),n,1);
            b = log(1-PI);
            B1 = log(PI) - b;
            %     B1(B1==-inf) = 0;
            int_iter = 1;
            diff = int_tol + 1;
            while (int_iter<int_maxIter) && (diff>int_tol)
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
                    diff = max(max(abs(Tau-old_Tau)));
                end
                if int_iter == 1
                    diff = max(max(abs(Tau-old_Tau)));
                end
                int_iter = int_iter + 1;
            end
            
        case 'Poisson'
            tm1 = repmat(log(Alpha),n,1);
            int_iter = 1;
            diff = int_tol + 1;
            while (int_iter<int_maxIter) && (diff>int_tol)
                old_Tau = Tau;
                tm2 = X*Tau*log(PI');
                tm3 = Tau*PI';
                tm3 = tm3 - repmat(sum(tm3,1),n,1);
                tm4 = log(factorial(X))*repmat(sum(Tau,2),1,Q);
                switch NetType
                    case 'undirected'
                        Tau = tm1 + 0.5*(tm2 + tm3 - tm4);
                    case 'directed'
                        Tau = tm1 + (tm2 + tm3 - tm4);
                end
                Tau(Tau<log_min) = log_min;
                Tau = Tau-repmat(max(Tau,[],2),1,Q);
                Tau = exp(Tau);
                Tau = Tau./repmat(sum(Tau,2),1,Q);
                Tau(Tau<prune) = 0;
                if int_iter > 1
                    diff = max(max(abs(Tau-old_Tau)));
                end
                if int_iter == 1
                    diff = max(max(abs(Tau-old_Tau)));
                end
                int_iter = int_iter + 1;
            end
    end
    out_iter = out_iter + 1;
end

time = toc;

%Est_Label
if Q == 1
    Est_Glabel = ones(1,n);
else
    [~,Est_Glabel] = max(Tau,[],2);
end