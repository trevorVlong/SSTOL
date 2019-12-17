%% Preamble
% Trevor Long
% 7 Sept, 2019
% Polynomial fit (looks like least squares regression)

function [avec,res,xmod1,xmod2,ymod,Yi] = polyregression(xdata1,xdata2,ydata,odr)

%% column check
sz1 = size(xdata1);
sz2 = size(xdata2);
sz3 = size(ydata);

if sz1(2) > 1 && sz1(1) == 1
    xdata1 = xdata1';
elseif sz1(2)> 1 && sz1(1) > 1
    error('input is not a vector')
end    
if sz2(2) > 1 && sz2(1) == 1
    xdata2 = xdata2';
elseif sz2(2)> 1 && sz2(1) > 1
    error('input is not a vector')
end    
if sz3(2) > 1 && sz3(1) == 1
    ydata = ydata';
elseif sz3(2)> 1 && sz3(1) > 1
    error('input is not a vector')
end    

%% correct order check
if odr == 1
    fprintf(1,'f(Dcj) = a1 * Dcj\n');
elseif odr == 2
    fprintf(1,'f(Dcj) = a0 + a1 * Dcj\n');
elseif odr == 3
    fprintf(1,'f(Dcj) = a0 + a1 * Dcj + a3*Dcj^2\n');
else
    error('incorrect order, 2nd or 3rd order only')
end

%% fitting process

xmod1 = linspace(0,20,21)';
xmod2 = linspace(-10,40,51)';
%outlier setup
% n_out = 1;
% N = 2;
% stdv = zeros
% iter = 0;
% while n_out ~= 0
    
    
    if odr == 1
        xmod1 = mean(xdata1);
        
        X = zeros(length(xdata1),6);
        Y = ydata;
        X(:,1) = 1;
        X(:,2) = xdata1;
        X(:,3) = xdata2;
        X(:,4) = xdata1.*xdata2;
        %X(:,5) = xdata2.^2;
        %X(:,6) = xdata1.*xdata2.^2;
        X(:,5) = xdata2.^3;
        X(:,6) = xdata1.*xdata2.^3;

        avec = (X'*X)\(X'*Y);
        
        Yi   = avec(1)+avec(2).*xdata1+avec(3).*xdata2+avec(4)*xdata1.*xdata2 ...
            ... %+ avec(5).*xdata2.^2 + avec(6)*xdata1.*xdata2.^2 ;
            + avec(5).*xdata2.^3 + avec(6)*xdata1.*xdata2.^3;
        
        for jj = 1: length(xmod1)
            ymod(:,jj) = avec(1) + avec(2).*xmod1(jj) + avec(3).*xmod2 + avec(4)*xmod1(jj).*xmod2 ...
                       ... + avec(5).*xmod2.^2 + avec(6)*xmod1(jj).*xmod2.^2 ;
                       + avec(5).*xmod2.^3 + avec(6)*xmod1.*xmod2.^3;
        end


    elseif odr == 2
        %2nd order best fit line
        X = zeros(length(xdata1),6);
        Y = ydata;
        X(:,1) = 1;
        X(:,2) = xdata1;
        X(:,3) = xdata2;
        X(:,4) = xdata1.*xdata2;
        X(:,5) = xdata2.^2;
        X(:,6) = xdata1.*xdata2.^2;

        avec = (X'*X)\(X'*Y);
        
        Yi   = avec(1)+avec(2).*xdata1+avec(3).*xdata2+avec(4)*xdata1.*xdata2 ...
            +avec(5).*xdata2.^2 + avec(6)*xdata1.*xdata2.^2;
        
        for jj = 1: length(xmod1)
            ymod(:,jj) = avec(1) + avec(2).*xmod1(jj) + avec(3).*xmod2 + avec(4)*xmod1(jj).*xmod2 ...
                       + avec(5).*xmod2.^2 + avec(6)*xmod1(jj).*xmod2.^2;
        end
        
    elseif odr == 3
        %2nd order best fit line
        X = zeros(length(xdata1),6);
        Y = ydata;
        X(:,1) = 1;
        X(:,2) = xdata1;
        X(:,3) = xdata2;
        X(:,4) = xdata1.*xdata2;
        X(:,5) = xdata2.^3;
        X(:,6) = xdata1.*xdata2.^3;

        avec = (X'*X)\(X'*Y);
        
        Yi   = avec(1)+avec(2).*xdata1+avec(3).*xdata2+avec(4)*xdata1.*xdata2 ...
            +avec(5).*xdata2.^3 + avec(6)*xdata1.*xdata2.^3;
        
        for jj = 1: length(xmod1)
            ymod(:,jj) = avec(1) + avec(2).*xmod1(jj) + avec(3).*xmod2 + avec(4)*xmod1(jj).*xmod2 ...
                       + avec(5).*xmod2.^3 + avec(6)*xmod1(jj).*xmod2.^3;
        end

    end
    
    
    res = (Yi-ydata); %residual for each point
    err2 = res.^2; %error squared
    
    
    
    
    %safety block
%     iter = iter+1;
%     if iter >= 10
%         fprintf('possible outliers\n')
%         break
%     end
% end

% build fitted domain




end
