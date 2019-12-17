function [CD,CL] = finitewing(cl,cd)
%converts 2D flight coefficients to 3D ones using standard corrections

%   wing basic parameters
    b  = 134; %inches
    c  = 18 ; %inches
    e  = .92; %non-dim spanwise efficiency
    AR = b/c;
    
    
    %current corrections, will change with funcitons later.
    CL  = .9*cl; %simple correction
    Cdi = CL.^2/(pi*AR*e); %assuming elliptic loading
    CD = cd+Cdi;
   
end

