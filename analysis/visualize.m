%plotbuild

% Trevor Long
% 13 Sept, 2019
close all;

%% choose run numbers to produce plots
% 
runmat = [46 53; ... 0  deg flaps
          54 60; ... 20 deg flaps
          61 68; ... 30 deg flaps
          70 77; ... 35 deg flaps
          81 88; ... 40 deg flaps
          89 96; ... 45 deg flaps
          103 109; ... 55 deg flaps
          97 102];   %60 deg flaps
%=========================================================================
%% build each subplot in loop
sz = size(runmat);

for plotn = 1:sz(1)
    [clmod,cxmod,cmmod,alfamod,Dcjmod,Dcjmodm,cli,cdi,cmi,alfamat,Dcjmat,Fdata,df] = rcoeff(runmat(plotn,1),runmat(plotn,2));
  % plot options
    ttl = sprintf('\x03b4_f = %02d\x00b0',df); %subplot labels
    al  = [-15,40];
    cll = [-2, 22];
    cxl = [-8, 4] ;
    cml = [-2,2];
    spm = 2;
    spn = 4;
  %build plot in subplot
  %cl-alpha
    visplot(1,plotn,2,4,al,cll,'\alpha','c_l' ,ttl,alfamat ,cli,Fdata);
    visplot(2,plotn,2,4,al,cxl,'\alpha','c_x' ,ttl,alfamat ,cdi,Fdata);
    visplot(3,plotn,2,4,al,cml,'\alpha','c_ml',ttl,alfamat,cmi,Fdata);
    visplot(4,plotn,2,4,cxl,cll,'c_x'   ,'c_l',ttl,cdi     ,cli,Fdata);
    
    fitplot(5,plotn,spm,spn,al,cll,'\alpha','c_l' ,ttl,alfamod,clmod,Dcjmod);
    fitplot(6,plotn,spm,spn,al,cxl,'\alpha','c_x' ,ttl,alfamod,cxmod,Dcjmod);
    %fitplot(7,plotn,spm,spn,al,cml,'\alpha','c_m' ,ttl,alfamod,cmmod,Dcjmod);
    fitplot(8,plotn,spm,spn,cxl,cll,'c_x','c_l' ,ttl,cxmod,clmod,Dcjmod);
    
end
F1 = figure(1);
F2 = figure(2);
F3 = figure(3);
F4 = figure(4);
F5 = figure(5);
F6 = figure(6);
F8 = figure(8);


%% 
F1.WindowState='maximized';
F2.WindowState='maximized';
F3.WindowState='maximized';
F4.WindowState='maximized';
F5.WindowState='maximized';
F6.WindowState='maximized';
F8.WindowState='maximized';

saveas(F1,'2D_cl_raw.eps','epsc')
saveas(F2,'2D_cx_raw.eps','epsc')
saveas(F3,'2D_cm_raw.eps','epsc')
saveas(F4,'2D_polar_raw.eps','epsc')
saveas(1,'2D_cl_raw.fig')
saveas(2,'2D_cx_raw.fig')
saveas(3,'2D_cm_raw.fig')
saveas(4,'2D_polar_raw.fig')

saveas(F5,'2D_cl_polyfit.eps','epsc')
saveas(F6,'2D_cx_polyfit.eps','epsc')
%saveas(7,'2D_cm_raw.epsc')
saveas(F8,'2D_polar_polyfit.eps','epsc')
saveas(F5,'2D_cl_polyfit.fig')
saveas(F6,'2D_cx_polyfit.fig')
%saveas(7,'2D_cm_raw.fig')
saveas(F8,'2D_polar_polyfit.fig')