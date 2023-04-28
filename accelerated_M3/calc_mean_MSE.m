function [MSE, Xt] = calc_mean_MSE(X)

[n, m] = size(X);
Xt = zeros(n, m);

count = 0;
for i = 1:m
    if (X(n, i) <= X(1, i))
        count = count + 1;
        Xt(:, count) = X(:, i);
    end;
end;

Xt = Xt(:, 1:count);
MSE = mean(Xt, 2);

return;