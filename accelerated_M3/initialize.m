function z = initialize(string,pos)
% Initialize the state vector with a realization of one of possible field types.
%% Inputs
% string - field type, string, the following values are supported:
%           'GB' - gaussian bumps
%           'Spike' - only one sensor holds 1, others hold 0
%           'IID' - i.i.d. Gaussian samples
%           'Slope' - sum of coordinates divided by 2
% pos - sensor coordinates
%% Outputs
% z - vector of initial values at the sensor nodes

n = size(pos, 1);

x = pos(:,1);
y = pos(:,2);

if strcmp(string,'GB')
    z = exp(-(8*(x-0.3)).^2 - (8*(y-0.4)).^2) ...
            + 7*exp(-(6.035*(x-0.65)).^2 - (6.035*(y-0.3)).^2) ...
            + 8*exp(-(10.065*(x-0.19)).^2 - (10.065*(y-0.19)).^2) ...
            + 25*exp(-(6.025*(x-0.5)).^2 - (6.025*(y-0.75)).^2);
elseif strcmp(string,'Spike')
    center = rand(2, 1);
    z = exp(-(50*(x-center(1))).^2 - (50*(y-center(2))).^2);
elseif strcmp(string,'IID')
    z=randn(n,1);
elseif strcmp(string,'Slope')
    z=(x+y)/2;
end

z = z - mean(z);
z = z / norm(z);

return;