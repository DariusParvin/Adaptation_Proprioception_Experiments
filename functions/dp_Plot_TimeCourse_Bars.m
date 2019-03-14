function dp_Plot_TimeCourse_Bars(bar_numbers, tableData, faceColor)
% First two columns of data are not for plotting

% bar_number = [1:4];


% Check that there is only one data point per subject
if  length(tableData.SN) ~= unique(tableData.SN)
    error('There is more than one datapoint per subj')
end

barVarNames = tableData.Properties.VariableNames;

for bi = 1:4;
    
    y1 = tableData.(barVarNames{bi+2});
    dpPlotBar2(bar_numbers(bi), y1, faceColor, 1);
    
end

end