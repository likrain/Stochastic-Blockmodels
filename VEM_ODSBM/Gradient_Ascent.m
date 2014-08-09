function W = Gradient_Ascent (W,A,Tau,Theta)

% Gradient Ascent Algorithm to Update W
% Author: Xiangyong Cao
% Date: 07/29/2014

% constant
Q = size(W,1);

maxIter = 100;
for q = 1:Q
    for l = 1:Q
        k = 1;
        % compute gradient at initial point
        dw = compute_gradient(W,A,Tau,q,l,Theta);
        % start Gradient Ascent iterations
        while(d'*d>=1e-10 && k<maxIter)
            % update w along the direction with a rate of 1/sqrt(k)
            w = w + 1/k*dw;
            
            % calculate gradient
            W(q,l) = w;
            dw = compute_gradient(W,A,Tau,q,l,Theta);
            
            % update k for next iteration
            k=k+1;
        end
    end
end
