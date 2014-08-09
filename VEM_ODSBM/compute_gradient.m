function dw = compute_gradient(W0,A,Tau,q,l,Theta)
% compute gradient at each iteration

% % dW: Q x Q
% tmp1 = Tau'*(A./(Tau*W*Tau'))*Tau;
% tmp2 = Tau'*(Theta*Theta')*Tau;
% dW = tmp1 - tmp2;

% dw: R
tm1 = Tau(:,q); tm2 = Tau(:,l);
tmp1 = tm1'*(A./(Tau*W*Tau'))*tm2;
tmp2 = tm1'*(Theta*Theta')*tm2;
dw = tmp1 - tmp2;
