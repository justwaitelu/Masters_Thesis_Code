clc
clear
close all

%grainsize caluclation
%addpath("EBSD\LucyWaite\2_3_2025")
%addpath("EBSD\LucyWaite\2_6_2025")
%addpath("EBSD\LucyWaite\2_7_2025")
%addpath("EBSD\LucyWaite\2_11_2025")

addpath("SEM\LucyWaite\3_12_2025")
addpath("SEM\LucyWaite\3_13_2025")
addpath("SEM\LucyWaite\3_14_2025")

%samples = ["A1","A2","A3","A4","A5","A6","A7","B1","B2","B3","B4","B5","B6","B7","C1","C2","C3","C4","C5","C6","C7","D1","D2","D3","D4","D5","D7","E1","E2","E3","E4","E5","F1","F2","F3","F4","F5"]';
samples = ["DUR", "D5R", "DHR", "DU8", "D58", "DH8", "FUR", "F5R", "FHR"];
averagegrainsize500 = zeros(length(samples),1); 
averagegrainsize5000 = zeros(length(samples),1); 

for i = 1:length(samples)
    samplefolder = samples(i);
    if exist("SEM\LucyWaite\3_12_2025\" + samplefolder) == 7, path = "SEM\LucyWaite\3_12_2025\" + samplefolder;
    elseif exist("SEM\LucyWaite\3_13_2025\" + samplefolder) == 7, path = "SEM\LucyWaite\3_13_2025\" + samplefolder;
    elseif exist("SEM\LucyWaite\3_14_2025\" + samplefolder) == 7, path = "SEM\LucyWaite\3_14_2025\" + samplefolder;
    else, disp("The sample folder could not be found")
    end

    addpath(path)

    grainsizetable500 = readtable("grainsize_500.txt");
    grainsizetable500 = grainsizetable500(1:25, :);

    % if strcmp(samples(i), "D5") | strcmp(samples(i), "D7")
    %     grainsizetable = readtable("Grain Size.txt");
    %     grainsizetable = grainsizetable(1:20, :);
    % end

    grainsizesmicron500 = table2array(grainsizetable500(:,1));
    percents500 = table2array(grainsizetable500(:,2));
    averagegrainsize500(i) = sum(grainsizesmicron500.*percents500); %microns

    if samplefolder ~= "FHR"
        grainsizetable5000 = readtable("grainsize_5000.txt");
        grainsizetable5000 = grainsizetable5000(1:25, :);
        grainsizesmicron5000 = table2array(grainsizetable5000(:,1));
        percents5000 = table2array(grainsizetable5000(:,2));
        averagegrainsize5000(i) = sum(grainsizesmicron5000.*percents5000); %microns
    end
   
    rmpath(path)

    % figure()
    % plot(grainsizesmicron, percents)
    % title(samples(i) + " Grain Size Distribution")
    % xlabel("Grain Size (micron)")
    % ylabel("Fraction")
end


compgrainsize = [samples', averagegrainsize500, averagegrainsize5000];
