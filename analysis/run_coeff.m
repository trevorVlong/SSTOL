%Trevor Long
%27 Aug 2019
% Wind Tunnel Rig data reduction
%%% REMEMBER TO AUTOMATE AoA SELECTION
clear ;
close all;


%==========================================================================
%% setup 
% open testlog
testlog = implog(); %testlog table

% import calibration
    %set up calibration matrix
    calfile = open("cal_0002.mat");
    calmat = cell2mat(struct2cell(calfile));
    

% import tare
    %import wind-off tare for test batch
    tarefile = open("tare_0017.mat");
    tarevec = cell2mat(struct2cell(tarefile))';
    %tarevec = [2.43 .22 3.09]';

% import raw data
    %import data (should be same number for each
    startnum = 97;
    endnum = 102;
    

    num = (endnum-startnum)+1;
    runs = {};
    skip = [];
    testlognum = [];
    for nn = 1:1:num
        % check point
        if testlog.GP(startnum+nn-1)== 'N' 
            skip = [skip (startnum+nn-1)];
            continue
        end
        
        file = sprintf("run_%04d.xlsx",nn+startnum-1);
        fileb = sprintf("run_%04d",nn+startnum-1);
        [a,b] = find(strcmp(testlog.filename(:),fileb) == 1);
        testlognum = [testlognum a];
        if isfile(file)
         % File exists.
            runs{nn} = importfile(file); 
        else
         %no file
        end
  
    end
    
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
%% Constants
r        = 287.05; % J/kgK specific ideal gas correction for air pv = rT rho = p/(rT)
airtemp  = 307; %kelvin temp in tunnel during run
q_t      = .22;       %torr tunnel dynamic pressure
static_p = 222; %static pressure on date
rho_pa   = 29.97 * 3386.39; %air density on day of tunnel run
rho      = 1.225;
A        = 0.139355; %wing area sq. meters
b        = 23.9375*0.0254; %meters
c        = 9*0.0254; %meters
Nu       = 1.81E-5; %SI Units
R        = 0.0635; %meters (prop radius)
rh       = 0.014; %meters  (hub radius)

% conversions
q        = 133.3224*q_t; %dynamic pressure in pascals
Vinf     = sqrt(2*q/rho); %freestream velocity m/s
%set zeros
vos      = calmat*tarevec; %V0 --offset/tare 0 values
sprintf('L = %04d \n D = %04d \n M = %04d',vos(1),vos(2),vos(3))


%==========================================================================
% pull out Voltages for forces and RPM in & out data to matrix

run1 = [];
for nn = 1:1:num

    if ismember(skip,(startnum+nn-1))
        continue
    end
    
    run1 = runs{nn};
    len = length(run1(:,1));
    for kk = 1:len
        cF_run(kk,:)  = (calmat*run1(kk,2:4)'-vos)/(q*A); %L,D,M
        RPM_in(kk,:)  = mean(run1(kk,9:10)); %RPM of motors 3
        PWM_out(kk,:) = run1(kk,8); %RPM
    end
    
    R_data(1:length(RPM_in),:,nn) = [cF_run, RPM_in,PWM_out]; %raw voltage and RPM data
    
end


%==========================================================================
%% reduction and analysis

PWM         = [1000,1050,1100,1150,1200,1300,1400,1500,1600];
lenp        = length(PWM);
Fdata   = [];
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
 odr = 2;
     n = numel(Fdata(:,1,:));
     alfamat = squeeze(Fdata(:,7,:));
     alfavec = reshape(alfamat,n,1); %get alfas into vector
     Dcjmat  = squeeze(Fdata(:,6,:));
     Dcjvec  = reshape(Dcjmat,n,1); %get Dcj """"
     clvec   = reshape(squeeze(Fdata(:,1,:)),n,1); %get cl """"
     cdvec   = reshape(squeeze(Fdata(:,2,:)),n,1); %get cd """""
     cmmat   = squeeze(Fdata(:,3,:)); %
     cmvec   = reshape(cmmat,n,1); %get cm """"""
%     
%     % input vector is(Dcj,alfa,ydata, fitting scheme)
     [cl_dcj(:,1),res,xmod1,xmod2,clymcj,cli]    = polyregression(Dcjvec,alfavec,clvec,3);
     [cd_dcj(:,1),res,Dcjmod,alfamod,cxymcj,cdi] = polyregression(Dcjvec,alfavec,cdvec,2);
     for n = 1:num
        [cm_dcj(:,n),res,Dcjmodcm(n),alfamodcm,cmymcj(:,n),cmi(:,n)] ...
         = polyregression(Dcjmat(n,:),alfamat(n,:),cmmat(n,:),1);
     end

     cli = reshape(cli,9,6);
     cdi = reshape(cdi,9,6);
     %cmi = reshape(cmi,9,8);


%==========================================================================
%% plotting
numplot = 16;
jmp = 4;

len = length(Fdata(1,1,:));
figure()
cmap = colormap('hot');

%Lift
for ii = 1:len
    DCJ = Fdata(:,6,ii);
    scatter(Fdata(:,7,ii),Fdata(:,1,ii),20,DCJ,'filled')
    hold on;
end
for n = 1:length(cli(:,1))
    plot(alfamat(n,:),cli(n,:))
    hold on
end

% for n = 1:jmp:numplot
% plot(alfamod,clymcj(:,n),'color',cmap(4*n,:));
% hold on
% end
%
colormap('hot')
cax = [0 15];
caxis(cax);
grid on
colorbar 
xlabel("\alpha");
ylabel("c_l");

%Drag
figure()


for ii = 1:len
    DCJ = Fdata(:,6,ii);
    scatter(Fdata(:,7,ii),Fdata(:,2,ii),20,DCJ,'filled')
    hold on    
end
for n = 1:length(cdi(:,1))
    plot(alfamat(n,:),cdi(n,:))
    hold on
end
% for n = 1:jmp:numplot
% plot(alfamod,cxymcj(:,n),'color',cmap(4*n,:));
% hold on
% end
len = length(Fdata(1,1,:));
colormap('hot')
caxis(cax);
grid on
colorbar 
xlabel("\alpha");
ylabel("c_x");

%Moment
figure()


for ii = 1:len
    DCJ = Fdata(:,6,ii);
    scatter(Fdata(:,7,ii),Fdata(:,3,ii),20,DCJ,'filled')
    hold on
    
end
% for n = 1:length(cmi(:,1))
%     plot(alfamat(n,:),cmi(:,n))
%     hold on
% end
for n = 1:num
plot(alfamodcm(:,1),cmymcj(:,n));
hold on
end
len = length(Fdata(1,1,:));
colormap('hot')
caxis(cax);
grid on
colorbar
xlabel('\alpha')
ylabel('c_m')

% Drag Polar
figure()


for ii = 1:len
    DCJ = Fdata(:,6,ii);
    scatter(Fdata(:,2,ii),Fdata(:,1,ii),20,DCJ,'filled')
    hold on;
end
for n = 1:length(cli(:,1))
    plot(cdi(n,:),cli(n,:))
    hold on
end

%  for n = 1:jmp:numplot
%  plot(cxymcj(:,n+1),clymcj(:,n),'color',cmap(4*n,:));
%  hold on
%  end
colormap('hot')
caxis(cax);
grid on
colorbar 
xlabel('c_x')
ylabel('c_l')