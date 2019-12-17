clear ;
close all;
startnum = 3;
endnum = 5;
%startnum = 6;
%endnum = 10;
%N = (endnum-startnum)+1;
iN = [10,11,12,13,14,15];
%iN = [6:11];
N = length(iN);
N_PWM = 6;
cL = zeros(N_PWM,N);
cL_J = zeros(N_PWM,N);
dCJ1 = zeros(N_PWM,N);
dCJ2 = zeros(N_PWM,N);
% VJ1 = zeros(N_PWM, N);
% VJ2 = zeros(N_PWM, N);
L  = zeros(N_PWM,N);
M  = zeros(N_PWM,N);
Vinf = zeros(1,N);
eta = zeros(N_PWM,N);
%[cL, cX, cM, dCJ1, dCJ2, L, X, M, PWM, AoA, Vinf,rho, VJ1, VJ2, dF, RPM1,RPM2,S] = reduce_test_point(filename)
for nn = 1:1:N
    run = sprintf("run_%04d",iN(nn)+100);
    if run == "run_0115"
        run = "run_0115_mod";
    end
    [cLs, ~, ~, dCJ1s, dCJ2s, Ls, ~, Ms, PWM_test, ~, Vinf(nn), rho, VJ1, VJ2,~,~,~,Sref] = reduce_test_point(run);
    dCJ = .5.*(dCJ1s(1:6)+dCJ2s(1:6));
    cL(:,nn) = cLs(1:6);
    VJ = (VJ1(1:6)+VJ2(1:6))*.5;
    L(:,nn) = Ls(1:6);
    M(:,nn) = Ms(1:6);
    cL_J(:,nn) = L(:,nn)./(.5.*rho.*transpose(VJ.^2).*Sref);
    
    PWM(:,nn) = PWM_test(1:6);
    eta(:,nn) = cL_J(:,nn)./dCJ';
end

figure()
hold on
for i = 1:N_PWM
    plot(Vinf, L(i,:))
end
xlabel('V_\infty (m/s)')
ylabel('L (N)')
title('Dimensional Lift')
legend(num2str(PWM(:,1)))


figure()
hold on
for i = 1:N_PWM
    plot(Vinf, M(i,:))
end
xlabel('V_\infty (m/s)')
ylabel('M (N*m)')
title('Dimensional Moment')
legend(num2str(PWM(:,1)))

figure()
hold on
for i = 1:N_PWM
    plot(Vinf, cL(i,:))
end
xlabel('V_\infty (m/s)')
ylabel('c_L')
title('C_L (relative to freestream)')
legend(num2str(PWM(:,1)))

figure()
hold on
for i = 1:N_PWM
    plot(Vinf, cL_J(i,:))
end
xlabel('V_\infty (m/s)')
ylabel('c_L')
title('C_L (relative to jet)')
legend(num2str(PWM(:,1)))

figure()
hold on
for i = 1:N_PWM
    plot(Vinf, eta(i,:))
end
xlabel('V_\infty (m/s)')
ylabel('c_L')

