function [cL, cX, cM, dCJ1, dCJ2, L, X, M, PWM, AoA, Vinf,rho, VJ1, VJ2, dF, RPM1,RPM2,S] = reduce_test_point(filename)
%REDUCE_TEST_POINT Returns the dimensional and dimensionless force values
%at a single test point, for each unique PWM setting. 
%   INPUTS
%   filename - name of run, associated with an entry in the testlog.xlsx
%               spreadsheet as well as a filename.xlsx data file.  For 
%               example, filename = 'run_0001'.

%   Outputs  - cL, 1xN_PWM vector of the 2D lift coeffficients. 
%              
%               PWM - vector of unique PWM values
%
%               N_PWM - length of PWM vector
%
%               
    if filename == "run_0115_mod"
        filename = "run_0115";
    end

    [q_tunnel, T_tunnel, p_amb,calfile, tarefile,AoA_hi, AoA_lo, dF] = get_testCondition(filename);
    tarevec = get_tarevec(tarefile);
    calfile = "cal_0002";
    calmat = cell2mat(struct2cell(open(strcat(calfile, ".mat"))));
    if filename == "run_0115"
        filename = "run_0115_mod";
    end
    run_data = importfile(strcat(filename, '.xlsx'));
    
    rawPWM = run_data(:,8);
    %The motors don't start reliably until a PWM signal of >1100 is commanded. 
    %There may be bad data where the PWM signal 1050 was commanded before
    %1100, so any values of 1050 that occur before 1100 will be thrown out.
    
    motorOff = 1;
    i_trim = [];
    j = 1;
    for i = 1:length(rawPWM)
        if motorOff
            if rawPWM(i) >= 1100
                motorOff = 0;
            else
                if rawPWM(i) > 1000 && rawPWM(i) < 1100
                    i_trim(j) = i;
                    j = j+1;
                end
            end
        end
    end
    
    if length(i_trim) > 0
        i_start = i_trim(1);
        i_end = i_trim(end);
        
        trimPWM = rawPWM([1:istart i_end:end]);
        trimdata = run_data([1:istart i_end:end],:);
    else
        trimPWM = rawPWM;
        trimdata = run_data;
    end
        
    % Constants (hardcoded geometry - change to configurable)
    R = 287.05; % J/kgK specific ideal gas correction for air pv = rT rho = p/(rT)
    b   =23.9375*0.0254;       %meters
    c   =9*0.0254;             %meters
    S   =b*c;                  %wing area sq. meters
    rp  =0.0635;    % meters
    rh  =0.0127;     % meters
    %Atmo corrections
    rho = p_amb./(R.*T_tunnel);
    Vinf = sqrt(2*q_tunnel./rho); %freestream velocity m/s
    
    
    
    [PWM, ~, iPWM] = unique(trimPWM);
    
    P = length(PWM);
    V_L = zeros(1,P);   %Lift load cell average voltage
    V_X = zeros(1,P);   %X load cell average voltage
    V_M = zeros(1,P);   %M load cell average voltage
    RPM1 = zeros(1,P);   %RPM1 average value
    RPM2 = zeros(1,P);   %RPM2 average value
    L = zeros(1,P);
    X = zeros(1,P);
    M = zeros(1,P);
    AoA = zeros(1,P);
    dCJ1 = zeros(1,P);
    dCJ2 = zeros(1,P);
    VJ1 = zeros(1,P);
    VJ2 = zeros(1,P);
    
    for p = 1:P
        BUFF = 2;       %Buffer to ignore any transients from motor spoolup or down. BUFF = 8 ~ corresponds to 1 second of data
        PWM_val = trimPWM(iPWM == p);
        subset = trimdata(iPWM == p,:);
        
        V_L(p) = mean(subset(BUFF:end-BUFF,2));
        V_X(p) = mean(subset(BUFF:end-BUFF,3));
        V_M(p) = mean(subset(BUFF:end-BUFF,4));
        
        RPM1(p) = mean(subset(BUFF:end-BUFF,9));
        RPM2(p) = mean(subset(BUFF:end-BUFF,10));
        l = 3/12*.3048;
        F = calmat*([V_L(p),V_X(p),V_M(p)]'-tarevec);
        L(p) = F(1);
        X(p) = F(2);
        M(p) = F(3)*l;
        
        if M(p) > 0
            AoA(p) = AoA_hi;
        else
            AoA(p) = AoA_lo;
        end
        if Vinf == 0
            CT1(p) = getrotorp(.01,RPM1(p));
            CT2(p) = getrotorp(.01,RPM2(p));
        else
            CT1(p) = getrotorp(Vinf,RPM1(p));
            CT2(p) = getrotorp(Vinf,RPM2(p));
        end
    end
    
    omega1 = RPM1.*2.*pi/60;
    omega2 = RPM2.*2.*pi/60;
    
    Adisk = pi*(rp^2-rh^2);
    T1 = CT1.*1/2.*rho.*(omega1*rp).^2*Adisk;
    T2 = CT2.*1/2.*rho.*(omega2*rp).^2*Adisk;
    
    VJ1=sqrt(Vinf^2+2*T1./(rho*Adisk));
    VJ2=sqrt(Vinf^2+2*T2./(rho*Adisk));
    Cpd=0;
    Nprop = 4;
    hd=Nprop*pi*(rp^2-rh^2)/b*sqrt(1-Cpd);
    
    CQ1=.5*(1+VJ1./Vinf)*hd/c;
    CQ2=.5*(1+VJ2./Vinf)*hd/c;
    CJ1=2.*CQ1.*VJ1./Vinf;
    CJ2=2.*CQ2.*VJ2./Vinf;
    
    dCJ1 = 2.*CQ1.*(VJ1./Vinf-Vinf./VJ1);
    dCJ2 = 2.*CQ2.*(VJ2./Vinf-Vinf./VJ2);
    
    cL = L./(q_tunnel*S);
    cX = X./(q_tunnel*S);
    cM = M./(q_tunnel*S);  

end

