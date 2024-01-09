function thresholds = thresholdfun(avals)
% calculates the threshold used in FPT analysis for n values of i between 
% amin and amax; returns vector of a values and corresponding thresholds

%Make global so our DE can change within for loops
global k;

% Storage of FPT threshold values
thresholds = zeros(length(avals),1);

%Calculate thresohld for each value of a
for k = 1:length(thresholds)
          
    %Starting values for Newtons Algorithm
    x0s = 1;
    
    %Make an empty array for storing maximum value of dx/dt; 
    zs = zeros(length(x0s),1);
    
    %Run Newtons algorithm initiatied at each x0
    for j = 1: length(zs)
        
        fun = @(x) DDE(x,avals); % make fun a function of x alone
        zs(j) = fzero(fun, x0s(j));
    end
    
    %Keep Only unique, positive roots
    zs = uniquetol(zs, .001);
    zs = zs(zs > 0);
    
    % Store the max value 
    thresholds(k) = max(zs);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function derDE = DDE(x,avals)
    %derivative of the deterministic ODE

    [r,K,h,q,~] = parameters();
    derDE = r-2.*r.*x/K-(avals(k).*q.*h.^q.*x.^(q-1))/((x.^q+h.^q).^2);
end

end

