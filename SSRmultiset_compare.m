clc
clear
close all

%filemanagement
filepaths = ["SSR/20241018/A","SSR/20241018/B","SSR/20241018/C", "SSR/20241018/D","SSR/20241018/F"];
filenames_A = ["RT tensile_100FCC","RT tensile_54FCC46HCP","RT tensile_100HCP","800C tensile_100FCC","800C tensile_54FCC46HCP","800C tensile_100HCP"];
filenames_C = ["RT tensile_100FCC","RT tensile_47FCC53HCP","RT tensile_100HCP","800C tensile_100FCC","800C tensile_47FCC53HCP","800C tensile_100HCP"];
filenames_D = ["RT tensile_100FCC","RT tensile_52FCC48HCP","RT tensile_100HCP","800C tensile_100FCC","800C tensile_52FCC48HCP","800C tensile_100HCP"];
filenames_B = ["RT tensile_100FCC_nearconv",'','',"800C tensile_100FCC_nearconv",'','']; %spaces left for 50% FCC and HCP
filenames_F = ["RT tensile_100FCC","RT tensile_53FCC47HCP","RT tensile_100HCP","800C tensile_100FCC","800C tensile_53FCC47HCP","800C tensile_100HCP"];


for i = 1:length(filepaths)
    addpath(filepaths(i))
end

colors = ["#0072BD","#D95319","#EDB120","#7E2F8E","#77AC30","#A2142F"];
%% Set A

figure('Name','Set A - 400W-1000mm/s-100 micron')
for i = 1:length(filenames_A)
    tensile_dataA(i) = readSSR(filepaths(1) + '/' + filenames_A(i));
    SetA_Data(i) = processSSR(tensile_dataA(i), char(filepaths(1)), filenames_A(i));
    if ismember(i,1:3)
        subplot(1,2,1)
        plot(SetA_Data(i).strain*100,SetA_Data(i).stress_Pa/1000000,'LineWidth', 4)
        xlabel('Strain (%)')
        ylabel('Stress (MPa)')
        ylim([0,1250])
        xlim([0,42])
        hold on
    end
    if ismember(i,4:6)
        subplot(1,2,2)
        plot(SetA_Data(i).strain*100,SetA_Data(i).stress_Pa/1000000,'LineWidth', 4)
        title(filenames_A(i), Interpreter="none")
        xlabel('Strain (%)')
        ylabel('Stress (MPa)')
        ylim([0,650])
        xlim([0,14])
        hold on
    end
end

subplot(1,2,1)
title("Room Temperature Tensile")
legend("100% FCC","54% FCC, 46% HCP", "100% HCP",'Location','southeast')
subplot(1,2,2)
title("800C"+char(176) +" Tensile")
legend("100% FCC","54% FCC, 46% HCP", "100% HCP",'Location','southeast')

%% Set C

figure('Name','Set C - 350W-800mm/s-100 micron')
for i = 1:length(filenames_C)
    tensile_dataC(i) = readSSR(filepaths(3) + '/' + filenames_C(i));
    SetC_Data(i) = processSSR(tensile_dataC(i), char(filepaths(3)), filenames_C(i));
    if ismember(i,1:3)
        subplot(1,2,1)
        plot(SetC_Data(i).strain*100,SetC_Data(i).stress_Pa/1000000,'LineWidth', 4)
        xlabel('Strain (%)')
        ylabel('Stress (MPa)')
        ylim([0,1250])
        xlim([0,42])
        hold on
    end
    if ismember(i,4:6)
        subplot(1,2,2)
        plot(SetC_Data(i).strain*100,SetC_Data(i).stress_Pa/1000000,'LineWidth', 4)
        title(filenames_C(i), Interpreter="none")
        xlabel('Strain (%)')
        ylabel('Stress (MPa)')
        ylim([0,650])
        xlim([0,14])
        hold on
    end
end

subplot(1,2,1)
title("Room Temperature Tensile")
legend("100% FCC","47% FCC, 53% HCP", "100% HCP",'Location','southeast')
subplot(1,2,2)
title("800C"+char(176) +" Tensile")
legend("100% FCC","47% FCC, 53% HCP", "100% HCP",'Location','southeast')

%% Set D
figure('Name','Set D - 350W-800mm/s-125 micron')
for i = 1:length(filenames_D)
    tensile_dataD(i) = readSSR(filepaths(4) + '/' + filenames_D(i));
    SetD_Data(i) = processSSR(tensile_dataD(i), char(filepaths(4)), filenames_D(i));
    if ismember(i,1:3)
        subplot(1,2,1)
        plot(SetD_Data(i).strain*100,SetD_Data(i).stress_Pa/1000000,'LineWidth', 4)
        xlabel('Strain (%)')
        ylabel('Stress (MPa)')
        ylim([0,1250])
        xlim([0,42])
        hold on
    end
    if ismember(i,4:6)
        subplot(1,2,2)
        plot(SetD_Data(i).strain*100,SetD_Data(i).stress_Pa/1000000,'LineWidth', 4)
        title(filenames_D(i), Interpreter="none")
        xlabel('Strain (%)')
        ylabel('Stress (MPa)')
        ylim([0,650])
        xlim([0,14])
        hold on
    end
end

subplot(1,2,1)
title("Room Temperature Tensile")
legend("100% FCC","52% FCC, 48% HCP", "100% HCP",'Location','southeast')
subplot(1,2,2)
title("800C"+char(176) +" Tensile")
legend("100% FCC","52% FCC, 48% HCP", "100% HCP",'Location','southeast')


%% Set B
figure('Name','NEAR CONV - Set B - 400W-800mm/s-100 micron')
for i = 1:length(filenames_B)
    if ismember(i,[1,4])
        tensile_dataB(i) = readSSR(filepaths(2) + '/' + filenames_B(i));
        SetB_Data(i) = processSSR(tensile_dataB(i), char(filepaths(2)), filenames_B(i));
    end

    %Bring back if all 6 data files are available
    % if ismember(i,1:3)
    %     subplot(1,2,1)
    %     plot(SetB_Data(i).strain*100,SetB_Data(i).stress_Pa/1000000,'LineWidth', 4)
    %     xlabel('Strain (%)')
    %     ylabel('Stress (MPa)')
    %     ylim([0,1250])
    %     xlim([0,42])
    %     hold on
    % end
    % if ismember(i,4:6)
    %     subplot(1,2,2)
    %     plot(SetB_Data(i).strain*100,SetB_Data(i).stress_Pa/1000000,'LineWidth', 4)
    %     title(filenames_B(i), Interpreter="none")
    %     xlabel('Strain (%)')
    %     ylabel('Stress (MPa)')
    %     ylim([0,650])
    %     xlim([0,14])
    %     hold on
    % end
end
plot(SetB_Data(1).strain*100,SetB_Data(1).stress_Pa/1000000,'LineWidth', 4)
xlabel('Strain (%)')
ylabel('Stress (MPa)')
ylim([0,1250])
xlim([0,42])
hold on
plot(SetB_Data(4).strain*100,SetB_Data(4).stress_Pa/1000000,'LineWidth', 4)
legend("Room Temp", "800C" ,'Location','southeast')
title('100% FCC')

% subplot(1,2,1)
% title("Room Temperature Tensile")
% legend("100% FCC","52% FCC, 48% HCP", "100% HCP",'Location','southeast')
% subplot(1,2,2)
% title("800C"+char(176) +" Tensile")
% legend("100% FCC","52% FCC, 48% HCP", "100% HCP",'Location','southeast')


%% Set F
figure('Name','Set F - 250W-400mm/s-100 micron')
for i = 1:length(filenames_F)
        tensile_dataF(i) = readSSR(filepaths(5) + '/' + filenames_F(i));
        SetF_Data(i) = processSSR(tensile_dataF(i), char(filepaths(5)), filenames_F(i));
    if ismember(i,1:3)
        subplot(1,2,1)
        plot(SetF_Data(i).strain*100,SetF_Data(i).stress_Pa/1000000,'LineWidth', 4)
        xlabel('Strain (%)')
        ylabel('Stress (MPa)')
        ylim([0,1250])
        xlim([0,42])
        hold on
    end
    if ismember(i,4:6)
        subplot(1,2,2)
        plot(SetF_Data(i).strain*100,SetF_Data(i).stress_Pa/1000000,'LineWidth', 4)
        title(filenames_F(i), Interpreter="none")
        xlabel('Strain (%)')
        ylabel('Stress (MPa)')
        ylim([0,650])
        xlim([0,14])
        hold on
    end
end

subplot(1,2,1)
title("Room Temperature Tensile")
legend("100% FCC","53% FCC, 47% HCP", "100% HCP",'Location','southeast')
subplot(1,2,2)
title("800C"+char(176) +" Tensile")
legend("100% FCC","53% FCC, 47% HCP", "100% HCP",'Location','southeast')

%% Percent plots
figure('Name', '100% FCC Tests')
subplot(1,2,1)
plot(SetA_Data(1).strain*100,SetA_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetB_Data(1).strain*100,SetB_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(2))
plot(SetC_Data(1).strain*100,SetC_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(1).strain*100,SetD_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(1).strain*100,SetF_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("Room Temperature Tensile")
legend("Set A", "NEAR CONV - Set B","Set C", 'Set D','Set F','Location','southeast')
ylim([0,1250])
xlim([0,42])

subplot(1,2,2)
plot(SetA_Data(4).strain*100,SetA_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetB_Data(4).strain*100,SetB_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(2))
plot(SetC_Data(4).strain*100,SetC_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(4).strain*100,SetD_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(4).strain*100,SetF_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("800C"+char(176) +" Tensile")
legend("Set A", "NEAR CONV - Set B","Set C", 'Set D','Set F','Location','southeast')
ylim([0,650])
xlim([0,14])

figure('Name', '~50% FCC Tests')
subplot(1,2,1)
plot(SetA_Data(2).strain*100,SetA_Data(2).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetC_Data(2).strain*100,SetC_Data(2).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(2).strain*100,SetD_Data(2).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(2).strain*100,SetF_Data(2).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("Room Temperature Tensile")
legend("Set A", "Set C", 'Set D','Set F','Location','southeast')
ylim([0,1250])
xlim([0,42])

subplot(1,2,2)
plot(SetA_Data(5).strain*100,SetA_Data(5).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetC_Data(5).strain*100,SetC_Data(5).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(5).strain*100,SetD_Data(5).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(5).strain*100,SetF_Data(5).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("800C"+char(176) +" Tensile")
legend("Set A", "Set C", 'Set D','Set F','Location','southeast')
ylim([0,650])
xlim([0,14])


figure('Name', '100% HCP Tests')
subplot(1,2,1)
plot(SetA_Data(3).strain*100,SetA_Data(3).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetC_Data(3).strain*100,SetC_Data(3).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(3).strain*100,SetD_Data(3).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(3).strain*100,SetF_Data(3).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("Room Temperature Tensile")
legend("Set A", "Set C", 'Set D','Location','southeast')
ylim([0,1250])
xlim([0,42])

subplot(1,2,2)
plot(SetA_Data(6).strain*100,SetA_Data(6).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetC_Data(6).strain*100,SetC_Data(6).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(6).strain*100,SetD_Data(6).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(6).strain*100,SetF_Data(6).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("800C"+char(176) +" Tensile")
legend("Set A", "Set C", 'Set D','Set F','Location','southeast')
ylim([0,650])
xlim([0,14])

%% Temperature Plots
figure('Name', 'Room Temperature Tests')
subplot(1,3,1)
plot(SetA_Data(1).strain*100,SetA_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetB_Data(1).strain*100,SetB_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(2))
plot(SetC_Data(1).strain*100,SetC_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(1).strain*100,SetD_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(1).strain*100,SetF_Data(1).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("100% FCC")
legend("Set A", "NEAR CONV - Set B","Set C", 'Set D','Set F','Location','southeast')
ylim([0,1250])
xlim([0,42])

subplot(1,3,2)
plot(SetA_Data(2).strain*100,SetA_Data(2).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetC_Data(2).strain*100,SetC_Data(2).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(2).strain*100,SetD_Data(2).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(2).strain*100,SetF_Data(2).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("~50%/50% FCC/HCP")
legend("Set A (54%FCC)", "Set C (47%FCC)", 'Set D (52%FCC)','Set F (53%FCC)','Location','southeast')
ylim([0,1250])
xlim([0,42])

subplot(1,3,3)
plot(SetA_Data(3).strain*100,SetA_Data(3).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetC_Data(3).strain*100,SetC_Data(3).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(3).strain*100,SetD_Data(3).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(3).strain*100,SetF_Data(3).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("100% HCP")
legend("Set A", "Set C", 'Set D','Set F','Location','southeast')
ylim([0,1250])
xlim([0,42])


figure('Name', "800C"+char(176) +" Tests")
subplot(1,3,1)
plot(SetA_Data(4).strain*100,SetA_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetB_Data(4).strain*100,SetB_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(2))
plot(SetC_Data(4).strain*100,SetC_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(4).strain*100,SetD_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(4).strain*100,SetF_Data(4).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("100% FCC")
legend("Set A", "NEAR CONV - Set B", "Set C", 'Set D','Set F','Location','southeast')
ylim([0,650])
xlim([0,14])

subplot(1,3,2)
plot(SetA_Data(5).strain*100,SetA_Data(5).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetC_Data(5).strain*100,SetC_Data(5).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(5).strain*100,SetD_Data(5).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(5).strain*100,SetF_Data(5).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("~50%/50% FCC/HCP")
legend("Set A (54%FCC)", "Set C (47%FCC)", 'Set D (52%FCC)','Set F (53%FCC)','Location','southeast')
ylim([0,650])
xlim([0,14])

subplot(1,3,3)
plot(SetA_Data(6).strain*100,SetA_Data(6).stress_Pa/1000000,'LineWidth', 4,'Color', colors(1))
hold on
plot(SetC_Data(6).strain*100,SetC_Data(6).stress_Pa/1000000,'LineWidth', 4,'Color', colors(3))
plot(SetD_Data(6).strain*100,SetD_Data(6).stress_Pa/1000000,'LineWidth', 4,'Color', colors(4))
plot(SetF_Data(6).strain*100,SetF_Data(6).stress_Pa/1000000,'LineWidth', 4,'Color', colors(5))
title("100% HCP")
legend("Set A", "Set C", 'Set D','Set F','Location','southeast')
ylim([0,650])
xlim([0,14])
