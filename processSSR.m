function output = processSSR(tensile_data, filepath, filename)
    datasettitle = "Set "+filepath(end)+" - "+filename;

    output.length_mm = tensile_data.length;
    strain = tensile_data.Displacement_mm_./output.length_mm;
    output.area_msq = pi*(tensile_data.width/2000)^2; %mm diameter to m radius
    stress = tensile_data.Force_n_/output.area_msq;

    output.ultimatestress_MPa = max(stress./1000000);
    output.fracturestrain = max(strain(tensile_data.Force_n_>500));

    stress = stress(strain<output.fracturestrain);
    strain = strain(strain<output.fracturestrain);

    output.stress_Pa = stress;
    output.strain = strain;
    
    % figure()
    % plot(tensile_data.elapsedTime, tensile_data.Displacement_mm_)
    % title(datasettitle + ' - Displacement vs Time')
    % xlabel('Elapsed Time (s)')
    % ylabel('Displacement (mm)')
    % 
    % figure()
    % plot(tensile_data.elapsedTime, tensile_data.Force_n_)
    % title(datasettitle + ' - Force vs Time')
    % xlabel('Elapsed Time (s)')
    % ylabel('Force (N)')
    % 
    % figure()
    % plot(tensile_data.Displacement_mm_,tensile_data.Force_n_)
    % title(datasettitle + ' - Force vs Displacement')
    % xlabel('Displacement (mm)')
    % ylabel('Force (N)')
    
    %Figure for stress and strain
    % figure()
    % plot(strain*100,stress/1000000)
    % title(datasettitle + ' - Stress vs Strain', Interpreter="none")
    % xlabel('Strain (%)')
    % ylabel('Stress (MPa)')
    % xlim([0,(output.fracturestrain+0.1)*100])
    % %ylim([0,600])
    % hold on
    
    straightelastic_region = find(strain < 0.01);
    elasticfit = polyfit(strain(straightelastic_region), stress(straightelastic_region),1);
    elasticstrainfit = 0:0.00001:0.03;
    elasticstressfit = polyval(elasticfit, elasticstrainfit);
    maxstress = max(stress);
    
    limitedelastic = find(elasticstressfit < maxstress);
    elasticstressfit = elasticstressfit(limitedelastic)';
    elasticstrainfit = elasticstrainfit(limitedelastic)';
    
    %plot((elasticstrainfit + 0.002*ones(size(elasticstressfit))).*100, elasticstressfit/1000000)

    [xi,yi] = polyxpoly(strain, stress, (elasticstrainfit + 0.002*ones(size(elasticstressfit))),elasticstressfit);
    
    %plot(xi.*100,yi./1000000,'*')

    output.yieldstress_MPa = yi./1000000;
    output.yieldstrain = xi;
    output.E_GPa = elasticfit(1)/1000000000;

    %text((output.fracturestrain)*50, output.ultimatestress_MPa/2+40, "Young's Modulus Calculated as "+num2str(output.E_GPa)+" GPa",'FontSize',7)
    %text((output.fracturestrain)*50, output.ultimatestress_MPa/2+10, "Yield Stress Calculated as "+num2str(output.yieldstress_MPa)+" MPa",'FontSize',7)
    %text((output.fracturestrain)*50, output.ultimatestress_MPa/2-20, "Ultimate Stress Calculated as "+num2str(output.ultimatestress_MPa)+" MPa",'FontSize',7)

    save(filepath+"/"+datasettitle+'.mat', "output")

    T = table(stress,strain);
    writetable(T, filepath+"/"+'stress&strain ' + datasettitle+'.txt')
end
