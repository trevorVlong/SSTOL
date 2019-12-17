%3-Axis load cell calibration
%
%Trevor Long 29 July, 2019

%% enter weights


V0=[2.76,0.130,3.285]';

%weights
WM = .200*9.81; %weight rested on wing
D1 = .314*9.81; %dragweight hangar
D2 = (.907 + D1/9.81)*9.81; %drag weight

%first step
VL1=[2.823 .160 3.230]';
VM1=[2.837 .168 4.345]';
VD1=[2.756 .700 3.230]';
%second step
VL2=[2.900 .300 3.182]';
VM2=[2.900 .276 5.384]';
VD2=[2.737 2.065  3.117]';

l=0.0762; %meters, moment arm length
B1=zeros(3,3);
B2=zeros(3,3);

%%
B1(:,1)=(VL1-V0)/(-WM);
B1(:,2)=(VD1-V0)/D1;
B1(:,3)=-((VM1-V0)/(-WM)-B1(:,1))/l;

B2(:,1)=(VL2-V0)/(-2*WM);
B2(:,2)=(VD2-V0)/D2;
B2(:,3)=-((VM2-V0)/(-2*WM)-B2(:,1))/l;

B=0.5*(B1+B2);
C=inv(B);
C1 = inv(B1)
C2 = inv(B2)