function Tave = calc_Tave(X, epsilon)

[n, m] = size(X);
Tave = zeros(m, 1);

X = 10*log10(X) < epsilon;

for i = 1:m
    
    temp = find(X(:, i), 1, 'first');
    if (~isempty(temp))
        Tave(i) = find(X(:, i), 1, 'first') -1;
    else
        Tave(i) = n;
    end;
end;

return;