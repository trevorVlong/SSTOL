%Qprop comparison

%Static Values
S_vel = 0; %m/s
S_thr  = 17.85; %Newtons (stand in value)

%fli
F_vel = 6.2; %m/s
hd = 0.1778;
c  = 0.4572;
Vinf = F_vel;


%% filenames

m_file = sprintf('Scorpion'); %name of motor file
p_file = sprintf('an7x6-5');  %name of prop file
                                                                                        
 % qprop propfile motorfile 4.0      0      0        0.0        0          0.03     0   0
 % run   ' ' ''             V(m/s)  RPM     Volt    dBeta(deg) Thrust(N)   Torque Amp   Pele
 % 0 = Unspecified (dBeta means 0)
 % reads first one
 
 fpath = sprintf('../../../Desktop/Terminal\\ Programs/Qprop/bin');
 apath = sprintf('../../../Desktop/Terminal Programs/Qprop/bin');
 Mpath = sprintf('/Users/trevorlong/Dropbox (MIT)/Wind Tunnel Testing/Testing_Summer2019');
 fout1 = sprintf('out1.dat');
 fout2 = sprintf('out2.dat');
 
 addpath(apath);
 addpath(Mpath);
cmd1 = sprintf('cd %s; ./qprop %s %s %04d %04d %04d %04d %04d  >%s',fpath,p_file,m_file,S_vel,0,0, 0,S_thr,fout1);
status1 = system(cmd1);
[Vms_s,rpm_s,Dbeta_s,TN_s,QNm_s,PshaftW_s,Volts_s,Amps_s,effmot_s,effprop_s,adv_s,CT_s,CP_s,DVms_s,eff_s,Pelec_s,Pprop_s,cl_avg_s,cd_avg_s] = QPsread(fout1);



cmd2 = sprintf('cd %s; ./qprop %s %s %04d %04d %04d %04d %04d  >%s',fpath,p_file,m_file,F_vel,0,.57*Volts_s, 0,0,fout2);
[status2,out] = system(cmd2);
status2
[Vms,rpm,Dbeta,TN,QNm,PshaftW,Volts,Amps,effmot,effprop,adv,CT,CP,DVms,eff,Pelec,Pprop,cl_avg,cd_avg] = QPsread(fout2);
tbl = QPsreadp(fout2);
Wa_avg = 1/(tbl.radius(end)-tbl.radius(1))*trapz(tbl.radius(:),tbl.Wam(:));

Vj = 2*Wa_avg-F_vel;

Dcj = hd/c*(Vj^2/Vinf^2)*(Vinf/Vj+1)


