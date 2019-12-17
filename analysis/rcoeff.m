%Trevor Long
%11 Sept 2019
% Wind Tunnel Rig data reduction

function [clmod,cdmod,cmmod,alfamod,Dcjmod,Dcjmodm,cli,cdi,cmi,alfamat,Dcjmat,Fdata,df] = rcoeff(startnum,endnum)
% Hard-Coded values

%==========================================================================
    %% setup 
    % open testlog
    testlog = implog(); %testlog table
    

    % import raw data
        %import data (should be same number for each

        num = (endnum-startnum)+1;
        runs = {};
        skip1 = [];
        skip2 = [];
        testlognum = [];
        for nn = 1:1:num
            fileb = sprintf("run_%04d",nn+startnum-1);
            [a,b] = find(strcmp(testlog.filename(:),fileb) == 1);
            testlognum = [testlognum a];
            % check point
            if testlog.GP(startnum+nn-1)== 'N'
                skip1 = [skip1 (startnum+nn-1)];
                skip2 = [skip2 nn];
                continue
            end

            file = sprintf("run_%04d.xlsx",nn+startnum-1);


            if isfile(file)
             % File exists.
                runs{nn} = importfile(file);

            else
             %no file
            end

        end

        df      = testlog.FA(testlognum(1)); %flap angle
    % import calibration
        %set up calibration matrix
        calnum   = testlog.Cfile(startnum);
        calfile  = open(sprintf("%s.mat",calnum));
        calmat   = cell2mat(struct2cell(calfile));


    % import tare
        %import wind-off tare for test batch
        tarenum  = testlog.Tfile(startnum);
        tarefile = open(sprintf("%s.mat",tarenum));
        tarevec  = cell2mat(struct2cell(tarefile))';



%==========================================================================
    %columns description
    % 1: Time
    % 2: Lift
    % 3: Drag
    % 4: Moment
    % 5: Accel x
    % 6: Accel y
    % 7: Accel z
    % 8: pwm out
    % 9: rpm1
    % 10: rpm2
    % 11: ?
    % 12: ?
%==========================================================================
%test conditions
q   = testlog.DP(testlognum(1));
p   = testlog.AP(testlognum(1));
T   = mean([testlog.TT(testlognum(1)), testlog.TT(testlognum(end))]);

%==========================================================================
    %% Constants
    r        = 287.05; % J/kgK specific ideal gas correction for air pv = rT rho = p/(rT)
    A        = 0.139355; %wing area sq. meters
    b        = 23.9375*0.0254; %meters
    c        = 9*0.0254; %meters
    Nu       = 1.81E-5; %SI Units
    R        = 0.0635; %meters (prop radius)
    rh       = 0.014; %meters  (hub radius)

    % conversions
    rho      = p/(r*T);
    Vinf     = sqrt(2*q/rho); %freestream velocity m/s

    %set zeros
    vos      = calmat*tarevec; %V0 --offset/tare 0 values
    %sprintf('Tare\n Loffset = %04d \n Doffset = %04d \n Moffset = %04d\n',vos(1),vos(2),vos(3))


%==========================================================================
    % pull out Voltages for forces and RPM in & out data to matrix

    for nn = 1:1:num

        if ismember(skip1,(startnum+nn-1))
            continue
        end

        run1 = runs{nn};
        len = length(run1(:,1));
        for kk = 1:len
            cF_run(kk,:)  = (calmat*run1(kk,2:4)'-vos); %L,D,M
            RPM_in(kk,:)  = mean(run1(kk,9:10)); %RPM of motors 3
            PWM_out(kk,:) = run1(kk,8); %RPM
        end
        cF_run(:,1:2)                 = cF_run(:,1:2)/(q*c);
        cF_run(:,3)                 = cF_run(:,3)/(q*c^2);
        R_data(1:length(RPM_in),:,nn) = [cF_run, RPM_in,PWM_out]; %raw voltage and RPM data

    end


%==========================================================================
    %% reduction and analysis

    PWM         = [1000,1050,1100,1150,1200,1300,1400,1500,1600];
    lenp        = length(PWM);
    Fdata       = [];
    for nn = 1:1:num
        for pp = 1:lenp
            [I,a] = find(R_data(:,:,nn) == PWM(pp));
            if isempty(I)
                continue
            else
                for kk = I
                    %assign first 5 columns (Lift, Drag, Moment, RPM, PWM)
                    Fdata(pp,1:5,nn) = [mean(R_data(I,1,nn)),mean(R_data(I,2,nn)),mean(R_data(I,3,nn)),mean(R_data(I,4,nn)),mean(R_data(I,5,nn))];
                    RPM                  = Fdata(pp,4,nn); %assign RPM
                    omega                = RPM*2*pi/60; %get angular velocity
                    CT                   = getrotorp(Vinf,RPM); %get thrust coefficient via qprop data
                    T                    = CT*1/2*rho*(omega*R)^2*pi*R^2; %thrust from motor
                    Vj                   = sqrt(Vinf^2+2*T/(rho*pi*R^2)); %jet velocity behin motor
                    Cpd                  = 0; %?
                    hd                   = c*pi*((R/c)^2-(rh/c)^2)*4/(b/c)*sqrt(1-Cpd);%velocity disk hight
                    CQ                   = 1/2*(1+Vj/Vinf)*hd/c; %?
                    CJ                   = 2*CQ*Vj/Vinf; %jet velocity coefficient

                    %handle omega exceptions
                    if omega>0
                        delta_CJ         = 2*CQ*(Vj/Vinf-Vinf/Vj); %jet velocity excess coefficient
                    elseif isnan(omega)
                        delta_CJ         = NaN; %jet velocity excess coefficient
                    else
                        delta_CJ         = 0; %jet velocity excess coefficient
                    end


                    Fdata(pp,6,nn)   = delta_CJ;
                end
            end
        end
    end

%==========================================================================
    %% get AoA
    AoAhi = testlog.AoAhi(testlognum);
    AoAlo = testlog.AoAlo(testlognum);

    for kk = 1:num
        
        if ismember(kk,skip2)
            continue
        end

        for jj = 1:length(Fdata(:,1,kk))
           if  Fdata(jj,3,kk) < 0
               Fdata(jj,7,kk) = AoAlo(kk);
           else
               Fdata(jj,7,kk) = AoAhi(kk);
           end
        end
        Fdata(:,:,kk)

    end


%==========================================================================
    %% curve fitting

    % %fitting coeffs to DCJ
         n       = numel(Fdata(:,1,:));
         alfamat = squeeze(Fdata(:,7,:));
         alfavec = reshape(alfamat,n,1); %get alfas into vector
         Dcjmat  = squeeze(Fdata(:,6,:));
         Dcjvec  = reshape(Dcjmat,n,1); %get Dcj """"
         clvec   = reshape(squeeze(Fdata(:,1,:)),n,1); %get cl """"
         cdvec   = reshape(squeeze(Fdata(:,2,:)),n,1); %get cd """""
         cmmat   = squeeze(Fdata(:,3,:)); %
         cmvec   = reshape(cmmat,n,1); %get cm """"""
         
         cl_raw = squeeze(Fdata(:,1,:));
         
    %     
    %     % input vector is(Dcj,alfa,ydata, fitting scheme)
         [cl_dcj(:,1),res,xmod1,xmod2,clymcj]    = polyregression(Dcjvec,alfavec,clvec,3);
         [cd_dcj(:,1),res,Dcjmod,alfamod,cxymcj] = polyregression(Dcjvec,alfavec,cdvec,2);
         for n = 1:length(Fdata(:,1,1))
            [cm_dcj(:,n),res,Dcjmodm(n),alfamodcm,cmymcj(:,n)] ...
             = polyregression(Dcjmat(n,:),alfamat(n,:),cmmat(n,:),1);
         end

         %cli     = reshape(cli,9,num);
         %cdi     = reshape(cdi,9,num);
         %cmi     = reshape(cmi,9,num);
         clmod   = clymcj;
         cdmod   = cxymcj;
         cmmod   = cmymcj;
         
         cli = squeeze(Fdata(:,1,:));
         cdi = squeeze(Fdata(:,2,:));
         cmi = squeeze(Fdata(:,3,:));
%==========================================================================
end