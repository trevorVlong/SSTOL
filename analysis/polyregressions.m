%% Preamble
% Trevor Long
% 7 Sept, 2019
% Polynomial fit (looks like least squares regression)

function [avec,res,xmod1,xmod2,ymod] = polyregressions(xdata1,xdata2,ydata,odr)

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
    fprintf(1,'f(Dcj) = a1 * Dcj');
elseif odr == 2
    fprintf(1,'f(Dcj) = a0 + a1 * Dcj');
elseif odr == 3
    fprintf(1,'f(Dcj) = a0 + a1 * Dcj + a3*Dcj^2');
else
    error('incorrect order, 2nd or 3rd order only')
end

%% fitting process

xmod1 = mean(xdata1);
xmod2 = linspace(-10,40,100)';
%outlier setup

    
    
    if odr == 1
        X = zeros(length(xdata1),6);
        Y = ydata;
        X(:,1) = 1;
        X(:,2) = xdata1;
        X(:,3) = xdata2;
        X(:,4) = xdata1.*xdata2;
        X(:,5) = xdata2.^2;
        X(:,6) = xdata1.*xdata2.^2;
        X(:,7) = xdata1.^3;
        X(:,8) = xdata2.^3;

        avec = (X'*X)/(X'*Y);

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

        avec = (X'*X)/(X'*Y);

    end
    
    %doin' some stats
    
    Yi   = avec(1)+avec(2)*xdata1+avec(3)*xdata2+avec(4)*xdata1.*xdata2+avec(5)*xdata2.^2+avec(6)*xdata1.*xdata2.^2 + ;
        ymod =  avec(1)+avec(2)*xmod1+avec(3)*xmod2+avec(4)*xmod1.*xmod2+avec(5)*xmod2.^2+avec(6)*xmod1.*xmod2.^2;
    sz = size(ymod);
    if sz(2) > 1
        ymod = ymod';
    end
    
% %     if odr == 3
% %         Yi = Yi + avec(4)*xdata.^3;
% %         ymod = ymod + avec(4)*x1mod.^3;
% %     end
    
    res = (Yi-ydata); %residual for each point
    err2 = res.^2; %error squared
    
    
    
    


% build fitted domain




end
