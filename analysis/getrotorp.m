function CTrpm = getrotorp(Vinf,RPM)
%get Thrust coefficient based on input V_infinity and motor rpm
    [V,rpm,Dbeta,T,Q,Pshaft,Volts,Amps,effmot,effprop,adv,ct,CP,DV,eff,Pelec,Pprop,clavg,cdavg] = QPimport('F40qp.dat');
    CTmat = [V,ct,rpm];  %thrust coefficient aligned with V
    
    % filter by Velocity then RPM
    
    %getting closest velocity
    [a,b] = min(abs(CTmat(:,1)-Vinf));
    val = CTmat(b,1);
    [a,b] = find( (CTmat(:,1) == val));
    rvec = CTmat(a,:); %reduced vector with only needed velocities
    
    %finding closest RPM
    [x,y] = min(abs(rvec(:,3)-RPM)); %find closest rpm
    CTrpm = rvec(y,2);

end