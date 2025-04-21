clear
close all

%Ni200 Data for quartz correction
%addpath("..\..\2023-2024\Spring\Research MATLAB Files\CoCrMO\Ni200\")
addpath("Ni200 Calibration\")

%Print Parameters of sample - change for each new plot set
laser_power = 350;
scan_speed = 800;
hatch_spacing = 125;
powderthickness = 30;
printset_fileend = strcat('C_',num2str(laser_power),'_',num2str(scan_speed),'_',num2str(hatch_spacing));

pathstr = "AM CoCrMo " + num2str(laser_power)+'_'+num2str(scan_speed)+'_'+num2str(hatch_spacing);
addpath(pathstr)
disp('Plotting Graphs and Data for ' + pathstr)
disp('Volumetric Energy Density (VED) = ')
power_density = laser_power/(hatch_spacing*scan_speed*powderthickness);
disp(power_density)

dimentionsdatasheet = readtable(pathstr+'\Phase Transition Data Sheets('+ num2str(laser_power)+'-'+num2str(scan_speed)+'-'+num2str(hatch_spacing)+').csv');

temparray = [780,800,810,820,830,840,850,860,870,880];
%holdfittemps = [700,750,780,800,810,820,830,840,850,860,870,880];

ten_times = zeros(1,length(temparray));
fifty_times = zeros(1,length(temparray));
ninty_times = zeros(1,length(temparray));
coeffs = zeros(length(temparray),2);
coeffs2 = zeros(length(temparray),1);
tenfit = zeros(1,length(temparray));
fiftyfit = zeros(1,length(temparray));
nintyfit = zeros(1,length(temparray));
gofs = struct('sse',{},'rsquare',{},'dfe',{},'adjrsquare',{},'rmse', {});
gofs2 = struct('sse',{},'rsquare',{},'dfe',{},'adjrsquare',{},'rmse', {});
staticn = [];

for i = 1:length(temparray)
    filename = strcat(num2str(temparray(i)),printset_fileend);
    if isfile(pathstr+"\"+filename+"_T2.txt")
       datatable = readtable(pathstr+"\"+filename+"_T2");
    else
       datatable = readtable(pathstr+"\"+filename);
    end

    if isfile(strcat("Ni200 Calibration\Ni200_", num2str(temparray(i)),"C.txt"))
        Ni200_data = readtable(strcat("Ni200 Calibration\Ni200_", num2str(temparray(i)),"C.txt"));
    else, Ni200_data = [];
    end
    
    holdtempidx_dims = find(dimentionsdatasheet.Var3(:) == temparray(i));
    sample_length = mean(dimentionsdatasheet.Var5(holdtempidx_dims)./1000); %convert to meters

    [ten_times(i),fifty_times(i),ninty_times(i),tenfit(i),fiftyfit(i),nintyfit(i), coeffs(i,:),gofs(i), coeffs2(i,:), gofs2(i), staticn] = findtimes(datatable,Ni200_data, temparray(i), sample_length, laser_power, scan_speed, hatch_spacing);
end


figure()
%semilogx(ten_times, temparray,'*','MarkerSize',10)
%hold on
%ylim([750,900])
%semilogx(fifty_times, temparray,'*','MarkerSize',10)
%semilogx(ninty_times,temparray,'*','MarkerSize',10)

semilogx(tenfit./60, temparray,'*','MarkerSize',10)
hold on
ylim([750,900])
semilogx(fiftyfit./60, temparray,'*','MarkerSize',10)
semilogx(nintyfit./60,temparray,'*','MarkerSize',10)

% tenspline = fit((tenfit./60)', (temparray)', 'smoothingspline', 'SmoothingParam',0.07);
% plot(tenspline, tenfit./60, temparray)
% 
% fiftyspline = fit((fiftyfit./60)', (temparray)', 'smoothingspline', 'SmoothingParam',0.07);
% plot(fiftyspline, fiftyfit./60, temparray)
% 
% nintyspline = fit((nintyfit./60)', (temparray)', 'smoothingspline', 'SmoothingParam',0.07);
% plot(nintyspline, nintyfit./60, temparray)

title(num2str(laser_power) + "W, " + num2str(scan_speed)+ "mm/s " + num2str(hatch_spacing) + "µm, Transition Times vs Temps")
%[~, idx] = min(ninty_times);
ylabel(strcat('Aging Temperature (C',char(176),')'))
xlabel('Aging Time (min)')
xlim([10,1200])

[min90, idxmin90] = min(nintyfit);

save(pathstr+"/TTT_Output_Times.mat",'ten_times', 'fifty_times', 'ninty_times', 'temparray')
save(pathstr+"/k_n_fitcoeffs.mat",'coeffs')

%legend('Apx 10% HCP Times', 'Apx 50% HCP Times', 'Apx 90% HCP Times','Avrami Fit 10% Times','Avrami Fit 50% Times','Avrami Fit 90% Times', 'Location','southwest')
legend('Fitted 10% Times','Fitted 50% Times','Fitted 90% Times', 'Location','southwest')

figure()
%plot(temparray, coeffs(:,1))
idx = 9;
temparrayK = temparray + 273.*ones(1, length(temparray));
invtemparrayK = 1./temparrayK;
hold on
lnK = log(coeffs(:,1));
plot(invtemparrayK,lnK)
model1 = fit(invtemparrayK(1:idx)', lnK(1:idx), 'C - K*x','Start', [-10,5]);
%model2 = fit(invtemparrayK(idx+1:end)', lnK(idx+1:end), 'C - K/x','Start', [-10,5]);
hold on
plot(model1)
%plot(model2)
title("ln K vs Temperature")
xlabel("1/Temperature (1/K) ")
ylabel("ln(K)")

%dynamic vals
R = 8.31446261815324;
dynamicQ = model1.K/(R)
dynamiclnK = model1.C

figure()
plot(temparrayK, coeffs(:,2), '*')
xlabel("Temperature(K)")
ylabel('n')

ave_n = mean(coeffs(1:8,2))
aveGOF = mean(cell2mat({gofs.rsquare}))

figure()
%plot(temparray, coeffs(:,1))
idx = 7;
temparrayK = temparray + 273.*ones(1, length(temparray));
invtemparrayK = 1./temparrayK; 
hold on
lnK = log(coeffs(:,1));
plot(invtemparrayK,lnK)
model3 = fit(invtemparrayK(1:idx)', lnK(1:idx), 'C - K*x','Start', [-10,5]);
%model2 = fit(invtemparrayK(idx+1:end)', lnK(idx+1:end), 'C - K/x','Start', [-10,5]);
hold on
plot(model3)
%plot(model2)
title("ln K vs Temperature - STATIC N")
xlabel("1/Temperature (1/K) ")
ylabel("ln(K)")

staticQ = model3.K/(R)
staticlnK = model3.C


function [ten_time, fifty_time, ninty_time, tenfit,fiftyfit,nintyfit, coeffs, gofs, coeffs2, gofs2, ave_n_val] = findtimes(table, Ni200_data, holdtemp, barlength, power, speed, hatch)
    CoCrMo_times = cell2mat({table.Var1});
    CoCrMotemps = cell2mat({table.Var2});
    CoCrModisp = cell2mat({table.Var3}); %microns
    length_micron = barlength*10^6;
    thermal_expansion_coeff = mean((CoCrModisp./length_micron)./CoCrMotemps);
    smootheddisps = smooth(CoCrModisp,50);

    if isempty(Ni200_data) == 0
        Ni200_times = cell2mat({Ni200_data.Var1});
        Ni200_disps = cell2mat({Ni200_data.Var3});
        if length(Ni200_disps)>= length(smootheddisps)
            Ni200_times = Ni200_times(1:length(CoCrMo_times));
            Ni200_disps = (Ni200_disps(1:length(CoCrModisp))) - 50*ones(length(CoCrModisp),1);
        else
            difference = length(smootheddisps) -length(Ni200_disps);
            stepsize = Ni200_times(end)-Ni200_times(end-1);
            Ni200_times = [Ni200_times; Ni200_times(end) + stepsize.*(1:difference)'];
            Ni200_disps = [Ni200_disps; Ni200_disps(end).*ones(difference, 1)];
        end
        Ni200_disps = smooth(Ni200_disps, 100);
    else, Ni200_disps = zeros(length(smootheddisps),1);
    end

    

    figure()
    semilogx(CoCrMo_times, CoCrModisp,'-m', LineWidth=3)
    hold on
    semilogx(CoCrMo_times, smootheddisps,'-g', LineWidth=3)
    xlim([5,Inf])
    title(strcat(num2str(holdtemp),'C Displacement vs Time'))
    xlabel('Time (min)')
    ylabel('Displacement (µm)')

    [maxval_sm, maxidx_sm] = max(smootheddisps);

    phasetransition_maxshift = [700,150;750,70;780,20;800,20;810,15;820,12;830,12;840,10;850,10;860,10;870,5;880,5];
    adjustment_idx = find(phasetransition_maxshift(:,1) == holdtemp);
    adjustment = phasetransition_maxshift(adjustment_idx,2);

    transition_disps = smootheddisps((maxidx_sm-adjustment):end);
    transition_times = CoCrMo_times((maxidx_sm-adjustment):end);
    Ni200_disps = Ni200_disps((maxidx_sm-adjustment):end);

    [minval_sm, minidx_sm] = min(transition_disps);

    phasetransition_minshift = [700,0;750,0;780,50;800,0;810,0;820,0;830,0;840,0;850,0;860,0;870,0;880,5];
    adjustment_idx2 = find(phasetransition_minshift(:,1) == holdtemp);
    adjustment2 = phasetransition_minshift(adjustment_idx2,2);

    transition_times = transition_times(1:minidx_sm-adjustment2);
    transition_disps = transition_disps(1:minidx_sm-adjustment2);
    Ni200_disps = Ni200_disps(1:minidx_sm-adjustment2);
    transition_disps = smooth(transition_disps,500);

    Ni200_disps = Ni200_disps - min(Ni200_disps)*ones(length(Ni200_disps),1);
    %plot(transition_times, Ni200_disps)
    % 
    % mdlNi200 = strcat("(b*x)^2 + ", num2str(min(transition_times)));
    % Ni200fit = fit(Ni200_disps, transition_times,mdlNi200,'Start', 5);
    % coeffsNi200 = coeffvalues(Ni200fit);
    % 
    % Ni200fitdisps =sqrt(transition_times - min(transition_times)*ones(length(Ni200_disps),1))./coeffsNi200(1);

    %semilogx(transition_times, transition_disps-Ni200_disps,'-b', LineWidth=2)
    %transition_disps = transition_disps-Ni200_disps;
    %semilogx(transition_times, Ni200fitdisps,'-', LineWidth=1)

    spread = maxval_sm - minval_sm;

    tenperdisp = round(maxval_sm - 0.1*spread,1);
    fiftyperdisp = round(maxval_sm - 0.5*spread,1);
    nintyperdisp = round(maxval_sm - 0.9*spread,1);

    ten_times = transition_times((round(transition_disps,1) == tenperdisp));
    fifty_times = transition_times((round(transition_disps,1) == fiftyperdisp));
    ninty_times = transition_times((round(transition_disps,1) == nintyperdisp));

    ten_time = mean(ten_times);
    fifty_time = mean(fifty_times);
    ninty_time = mean(ninty_times);
    
    plot(ten_time, tenperdisp,'*','MarkerSize',10)
    plot(fifty_time, fiftyperdisp,'*','MarkerSize',10)
    plot(ninty_time, nintyperdisp,'*','MarkerSize',10)

    if holdtemp == 800
        figure()
        semilogx(CoCrMo_times, smooth(smootheddisps,100),'o', MarkerSize=3)
        hold on
        semilogx(transition_times, transition_disps,'-', LineWidth=4, Color="#D95319")
        semilogx(transition_times,  Ni200_disps, '-.', LineWidth=2, Color="#4DBEEE")
        semilogx(transition_times,  transition_disps - (Ni200_disps - min(Ni200_disps)*ones(length(Ni200_disps),1)),'-.', LineWidth=2, Color="#7E2F8E")
        xlim([5,Inf])
        ylim([-5,250])
        title(strcat(num2str(holdtemp),'C Displacement vs Time'))
        xlabel('Time (min)')
        ylabel('Elongation (µm)')
        legend("Smoothed Elongation vs Time", "Data Section Used for Avrami Curve Fitting", "Smoothed Ni200 Expansion", "Quartz Expansion Effect Removed")
    end

    mindisps = min(transition_disps);
    mintimes = min(transition_times);
    maxdisps = max(transition_disps);
    maxtimes = max(transition_times);

    %%"a*(1-exp(-(k*(x-l)^n)))+p"
    a = maxdisps-mindisps;
    p = mindisps;
    l = maxtimes-mintimes;

    adjustedtimes = (transition_times-mintimes.*ones(length(transition_times),1)).*60;
    adjusteddisps = ones(length(transition_disps),1)-(transition_disps-mindisps.*ones(length(transition_disps),1))./a;
    
    fiftyperadj_time = max(adjustedtimes)/2;

    %Dynamic N
    mdleqn = "1-exp(-(k*x)^n)";
    [mdl,gofs] = fit(adjustedtimes, adjusteddisps, mdleqn,'Start', [round(1/fiftyperadj_time,6),2.5]);
    coeffs = coeffvalues(mdl);
    mdltimes = (0:0.5:(round(max(adjustedtimes))+5000))';
    mdldisps = ones(length(mdltimes),1)-exp(-(coeffs(1).*mdltimes).^(coeffs(2)));
    figure()
    plot(adjustedtimes, adjusteddisps.*100,'o')
    hold on
    plot(mdltimes,mdldisps.*100,'r', 'LineWidth',3)
    title((num2str(holdtemp))+"C %HCP vs Time")
    xlabel("Time (s)")
    ylabel("HCP %")

    tenfitidx = round(mean(find(round(mdldisps,3) ==0.1)));
    fiftyfitidx = round(mean(find(round(mdldisps,3) ==0.5)));
    nintyfitidx = round(mean(find(round(mdldisps,3) ==0.9)));

    tenfit = mdltimes(tenfitidx)+mintimes*60;
    fiftyfit = mdltimes(fiftyfitidx)+mintimes*60;
    nintyfit = mdltimes(nintyfitidx)+mintimes*60;
    plot(mdltimes(tenfitidx),mdldisps(tenfitidx)*100,'*', 'MarkerSize',10,'LineWidth',1.5)
    plot(mdltimes(fiftyfitidx),mdldisps(fiftyfitidx)*100,'*','MarkerSize',10, 'LineWidth',1.5)
    plot(mdltimes(nintyfitidx),mdldisps(nintyfitidx)*100,'*', 'MarkerSize',10,'LineWidth',1.5)
    legend("Original Data (0-100% of the length change)","Avrami Fit", "10%", '50%','90%')


    %Static N, in order of 250_400_100, 350_800_100, 350_800_125,%400_1000_100
    ave_n = [1.535612140028430,1.683220972797012,1.771552332352303,1.760367704938110];

    if power == 250&&speed == 400&&hatch==100,    ave_n_val = ave_n(1);
    elseif power == 350&&speed == 800&&hatch==100,    ave_n_val = ave_n(2);
    elseif power == 350&&speed == 800&&hatch==125,    ave_n_val = ave_n(3);
    elseif power == 400&&speed == 1000&&hatch==100,    ave_n_val = ave_n(4);
    end

    mdleqn2 = "1-exp(-(k*x)^"+num2str(ave_n_val)+")";
    [mdl2,gofs2] = fit(adjustedtimes, adjusteddisps, mdleqn2,'Start', [round(1/fiftyperadj_time,6)]);
    coeffs2 = coeffvalues(mdl2);
    mdldisps2 = ones(length(mdltimes),1)-exp(-(coeffs2(1).*mdltimes).^ave_n_val);
    figure()
    plot(adjustedtimes, adjusteddisps.*100,'o')
    hold on
    plot(mdltimes,mdldisps2.*100,'r', 'LineWidth',3)
    title((num2str(holdtemp))+"C %HCP vs Time")
    xlabel("Time (s)")
    ylabel("HCP %")
    legend("Original Data (0-100% of the length change)", "Avrami Fit, fixed n")

end


% function [Ni200holddisp] = getNi200(holdtemp)
%     holdtempK = holdtemp + 273; %convert to kelvin;
% 
%     switch holdtemp    
%         case 700, data = readtable("Ni200-SHT-700C.txt");
%         case 750, data = readtable("Ni200-SHT-750C-100H.txt");
%         case 780, data = readtable("Ni200-SHT-780C-100H.txt");
%         case 800, data = readtable("Ni200-SHT-800C-100H.txt");
%         case 810, data = readtable("Ni200-SHT-810C-100H.txt");
%         case 820, data = readtable("Ni200-SHT-820C-100H.txt");
%         case 830, data = readtable("Ni200-SHT-830C-100H.txt");
%         case 840, data = readtable("Ni200-SHT-840C.txt");
%         case 850, data = readtable("Ni200-SHT-850C-100H.txt");
%         case 860, data = readtable("Ni200-SHT-860C.txt");
%         case 870, data = readtable("Ni200-SHT-870C.txt");
%         case 880, data = readtable("Ni200-SHT-880C.txt");
%     end
% 
%     Ni200alphafit = polyfit([973,1073,1173], [15.8,16.2,16.5],1);
%     Ni200_alpha = polyval(Ni200alphafit,holdtempK);
% 
%     Ni200_times = data.Var1;
%     Ni200_temps = data.Var2;
%     Ni200disps = data.Var3;
%     Ni200cleaneddisps = smooth(Ni200disps, 'sgolay');
%     [maxval, maxidx] = max(Ni200cleaneddisps);
%     lengthval = (maxval./holdtempK)./Ni200_alpha;
% 
%     holdidx = find(Ni200_temps>(holdtemp-1));
%     Ni200holddisp = Ni200cleaneddisps(holdidx);
%     Ni200holdtemp = Ni200_temps(holdidx);
%     Ni200holdtime = Ni200_times(holdidx);
%     Ni200thermalstrain = (Ni200holdtemp+ 273*ones(length(Ni200holdtemp),1)).*Ni200_alpha.*lengthval;
%     quartzisodisp = (Ni200holddisp-Ni200thermalstrain)+5*ones(length(Ni200holddisp),1);
%     quartzisodisp = smooth(quartzisodisp, 1000);
% 
%     figure()
%     plot(Ni200holdtime./max(Ni200holdtime), quartzisodisp)
%     title('holdtime v disp')
% 
%     ft2 = fittype(@(a,b,c,d,x) a*log(x-c).^b + d);
%     [holdfitval, holdfitgof] = fit(Ni200holdtime./max(Ni200holdtime), quartzisodisp,ft2,'Start',[1, 0.1, -0.9, 0]);
%     plot(holdfitval)
% end