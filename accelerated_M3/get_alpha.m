function [alp_mh, l2_mhM3] = get_alpha(l2_mh, l2_mh_est, theta)
% Calculate the optimal value of alpha, according to Theorem 1
%% Inputs
% l2_mh the eigenvalue of the foundational matrix W.
% l2_mh_est the estimate of the eigenvalue of the foundational matrix W if
% we have perfect knowledge of this egenvalue, set l2_mh_est = l2_mh.
% theta the vector of predictor weights.
%% Outputs
% alp_mh the value of alpha calculated, according to Theorem 1. If we have
% l2_mh_est = l2_mh then this corresponds to the optimal alpha identified
% in Theorem 1
% l2_mhM3 the value of the spectral radius of the prediction operator

alp_mh = (-((theta(3)-1)*l2_mh_est.^2 +theta(2)*l2_mh_est +2*theta(1) ) - 2*sqrt(theta(1)*(theta(1) +l2_mh_est.*(theta(2) +(theta(3)-1).*l2_mh_est))) )...
    ./ ((l2_mh_est*(theta(3)-1) + theta(2)).^2);

l2_mh_a = (1-alp_mh+alp_mh*theta(3))*l2_mh + alp_mh*theta(2);
l2_mhM3 = abs(1/2*(l2_mh_a + sqrt(l2_mh_a.^2 + 4*alp_mh*theta(1))));

return;