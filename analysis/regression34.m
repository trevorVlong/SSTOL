%% Preamble
% Trevor Long
% 5 Aug, 2019
% Description: completes 3rd or 4th order least-squares regression on two
% columns of data and removes outliers greater than 2 std from the best-fit
% line

function [avec,res,xmod,ymod] = regression34(xdata,ydata,odr)

%% column check
sz1 = size(xdata);
sz2 = size(ydata);

if sz1(2) > 1 && sz1(1) == 1
    xdata = xdata';
elseif sz1(2)> 1 && sz1(1) > 1
    error('input is not a vector')
end    

if sz2(2) > 1 && sz2(1) == 1
    ydata = ydata';
elseif sz2(2)> 1 && sz2(1) > 1
    error('input is not a vector')
end    

%% correct order check
if odr == 2
    
elseif odr == 3
    
else
    error('incorrect order, 2nd or 3rd order only')
end

%% fitting process

xmod = linspace(min(xdata),max(xdata),100);

%outlier setup
n_out = 1;
N = 2;
stdv = zeros
iter = 0;
while n_out ~= 0
    
    
    if odr == 3
        avec = zeros(4,1);
        X = zeros(length(xdata),4);
        Y = ydata;
        X(:,1) = 1;
        X(:,2) = xdata;
        X(:,3) = xdata.^2;
        X(:,4) = xdata.^3;

        avec = inv(X'*X)*X'*Y;

    elseif odr == 2
        %2nd order best fit line
        avec = zeros(3,1);
        X = zeros(length(xdata),3);
        Y = ydata;
        X(:,1) = 1;
        X(:,2) = xdata;
        X(:,3) = xdata.^2;

        avec = inv(X'*X)*X'*Y;

    end
    
    %doin' some stats
    
    Yi = avec(1) + avec(2)*xdata + avec(3)*xdata.^2;
    
    ymod = avec(1) + avec(2)*xmod + avec(3)*xmod.^2;
    
    if odr == 3
        Yi = Yi + avec(4)*xdata.^3;
        ymod = ymod + avec(4)*xmod.^3;
    end
    
    res = (Yi-ydata); %residual for each point
    err2 = res.^2; %error squared
    
    
    
    
    %safety block
    iter = iter+1;
    if iter >= 10
        fprintf('possible outliers\n')
        break
    end
end

% build fitted domain




end
