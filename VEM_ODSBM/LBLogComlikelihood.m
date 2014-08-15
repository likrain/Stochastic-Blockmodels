L = LBLogComlikelihood(A,Alpha,Theta,W,Tau)

S1 = log(Theta)*(sum(A,2)+sum(A,1)');
S2 = sum(sum((Tau*W*Tau').*(Theta'*Theta)));
tmp = repmat(Alpha,n,1);
S3 = sum(sum(Tau.*log(tmp)+(1-Tau).*log(1-tmp)));
S4 = sum(sum(Tau.*log(Tau)+(1-Tau).*log(1-Tau)));

L = S1 - S2 + S3 - S4;