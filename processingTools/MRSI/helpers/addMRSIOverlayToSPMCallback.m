function spmGlobalVariable = addMRSIOverlayToSPMCallback(spmGlobalVariable, functionName)
    arguments
        spmGlobalVariable (1, 1) struct
        functionName(:, 1) char
    end
    spmCallbackList = spmGlobalVariable.callback;
    if(ischar(spmCallbackList))
        %create cell array and add overaly() callback onto it
        callbackList = spmCallbackList;
        spmCallbackList = cell(1, 2);
        spmCallbackList{1} = callbackList;
        spmCallbackList{2} = functionName;
    else
        %append the overlay callback to the last spot in the cell array
        callbackStackSize = size(spmCallbackList, 2);
        spmCallbackList{callbackStackSize + 1} = functionName;
    end
    spmGlobalVariable.callback = spmCallbackList;
end