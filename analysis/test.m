F = [];
alfa = 0:5:30
for jj = 1:length(Dcj)
for ii = 1:length(alfa)
F(jj) = (a*Dcj(jj))'*[1 alfa(ii) alfa(ii)^3]';
end
plot(alfa,F)
hold on;
end