function [metabolitePlotArray, lcModelYIndexRange, lcModelXIndexRange] = makeLcModelHeatMap(lcModelTable, selectedMetaboliteIndex)
    lcModelYIndexRange(1) = min(lcModelTable.Row);
    lcModelYIndexRange(2) = max(lcModelTable.Row);
    lcModelXIndexRange(1) = min(lcModelTable.Col);
    lcModelXIndexRange(2) = max(lcModelTable.Col);

    numXIndecies = lcModelXIndexRange(2) - lcModelXIndexRange(1) + 1;
    numYIndecies = lcModelYIndexRange(2) - lcModelYIndexRange(1) + 1;
    metabolitePlotArray = zeros(numXIndecies, numYIndecies);
    selectedMetaboliteIndex = selectedMetaboliteIndex + 1;
    if(selectedMetaboliteIndex ~= 0)
        for iRow = 1:height(lcModelTable)
            xIndex = lcModelTable.Col(iRow);
            xIndex = xIndex - lcModelXIndexRange(1) + 1;
            yIndex = lcModelTable.Row(iRow);
            yIndex = yIndex - lcModelYIndexRange(1) + 1;

            metaboliteConcentration = lcModelTable(iRow, selectedMetaboliteIndex);
            metabolitePlotArray(xIndex, yIndex) = metaboliteConcentration.Variables;
        end
    else
        metabolitePlotArray = [];
    end
    fprintf('Current Metabolite is %s\n', lcModelTable.Properties.VariableNames{selectedMetaboliteIndex});
end