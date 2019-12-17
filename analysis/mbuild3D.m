%3DMapBuild

% Trevor Long
% 7 Oct, 2019
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
cl3D = [];
cx3D = [];
cm3D = [];
alfaq = linspace(-5,30,36);
Dcjq  = linspace(0,15,16);
dfvec = [0 20 30 35, 40, 45,55,60];
for matn = 1:sz(1)
    [clmod,cxmod,cmmod,alfamod,Dcjmod,Dcjmodm,cli,cdi,cmi,alfamat,Dcjmat,Fdata,df] = rcoeff(runmat(matn,1),runmat(matn,2));
    
    cmi = cmi';
    %build cell array to read each datapoint from
    method = 'linear';
    %method = 'nearest';
    if matn == 1
       Fcx_minor =  scatteredInterpolant(alfamat(1,:)',Dcjmat(1,:)',cdi(1,:)',method,'linear');
    end
    Fcl = scatteredInterpolant(alfamat(:),Dcjmat(:),cli(:),method,'linear');
    Fcx = scatteredInterpolant(alfamat(:),Dcjmat(:),cdi(:),method,'linear');
    Fcm = scatteredInterpolant(alfamat(:),Dcjmat(:),cmi(:),method,'linear');
    
    cl3D(:,:,matn) = Fcl({alfaq,Dcjq});
    cx3D(:,:,matn) = Fcx({alfaq,Dcjq});
    cm3D(:,:,matn) = Fcm({alfaq,Dcjq});
end
%%
if strcmp(method, 'nearest')
    save('cl3D_near.mat','cl3D');
    save('cx3D_near.mat','cx3D');
    save('cm3D_near.mat','cm3D');
    %save('paramat.mat','paramat');
    save('alfa_coord_near.mat','alfaq');
    save('Dcj_coord_near.mat','Dcjq');
    save('df_coord_near.mat','dfvec');
else
    save('cl3D.mat','cl3D');
    save('cx3D.mat','cx3D');
    save('cm3D.mat','cm3D');
    %save('paramat.mat','paramat');
    save('alfa_coord.mat','alfaq');
    save('Dcj_coord.mat','Dcjq');
    save('df_coord.mat','dfvec');
end