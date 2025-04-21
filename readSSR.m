function output = readSSR(filepath)
    datainfo = readtable(filepath + ".txt");
    datatable = readtable(filepath+'.dat');
    gettablenames = datatable.Properties.VariableNames;
    output = struct();
    for i = 1:length(gettablenames)
        output.(gettablenames{i}) = table2array(datatable(:,i));
    end
    minval = min(output.TimeStamp);
    elapsedTime = (output.TimeStamp - minval*ones(size(output.TimeStamp)));
    output.elapsedTime = elapsedTime; %in seconds

    lengthrow = find(ismember(table2array(datainfo(:,1)),'Length'))
    widthrow = find(ismember(table2array(datainfo(:,1)), 'Width'))
    output.length = table2array(datainfo(lengthrow,2))
    output.width = table2array(datainfo(widthrow,2))
end