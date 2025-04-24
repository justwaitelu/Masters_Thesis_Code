clc
clear
close all

power = [400 350 350 250];
speed = [1000 800 800 400];
hatch = [100 100 125 100];
power_density = power./(hatch.*speed.*30).*(1000000);
pathstr = [""; "";"";""];
temparray = [780,800,810,820,830,840,850,860,870,880];

for i = 1:length(power)
    pathstr(i) = "AM CoCrMo " + num2str(power(i))+'_'+num2str(speed(i))+'_'+num2str(hatch(i));
    addpath(pathstr(i))
    k_and_n_array(i,:) = load(strcat(pathstr(i),"\k_n_fitcoeffs.mat"));
end

figure()
plot(temparray, k_and_n_array(1).coeffs(:,1),"-*")
hold on
plot(temparray, k_and_n_array(2).coeffs(:,1),"-o")
plot(temparray, k_and_n_array(3).coeffs(:,1),"-+")
plot(temparray, k_and_n_array(4).coeffs(:,1),"-s")
legend("Set A","Set C", "Set D", "Set F")
xlabel("Temperature (C^\circ)")
ylabel("k-values")
xlim([770, 890]);
title("k vs Aging Temperature")

figure()
plot(temparray, k_and_n_array(1).coeffs(:,2),"*", 'MarkerSize',9)
hold on
plot(temparray, k_and_n_array(2).coeffs(:,2),"o", 'MarkerSize',9)
plot(temparray, k_and_n_array(3).coeffs(:,2),"+", 'MarkerSize',9)
plot(temparray, k_and_n_array(4).coeffs(:,2),"s", 'MarkerSize',9)
legend("Set A","Set C", "Set D", "Set F")
xlim([770, 890]);
xlabel("Temperature (C^\circ)")
ylabel("n-values")
title("n vs Aging Temperature")

tempK = temparray + 273*ones(1, length(temparray));
R = 8.31446261815324;
limit1 = 1:7;
limit2 = 1:8;

lnk = zeros(length(power_density), length(tempK));
figure()
hold on
for i = 1:length(power_density)
    if ismember(i, [1,3,4]),lim = limit1;
    elseif ismember(i, [2]), lim = limit2;  end
    invK = 1./(tempK(lim));
    lnk(i,lim) = log(k_and_n_array(i).coeffs(lim,1)');
    fitk(i,:) = polyfit(invK, lnk(i,lim), 1);
    k0(i) = exp(fitk(i,2));
    Q_R(i) = fitk(i,1);
    Q(i) = -Q_R(i).*R;
    plot(invK, lnk(i,lim))

    secondrange = max(lim)+1:length(tempK);
    lnk2 = log(k_and_n_array(i).coeffs(secondrange,1)');
    invK2 = 1./(tempK(secondrange));
    Q2 = 
    plot(invK2, lnk2)
end

