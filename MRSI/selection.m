%selection.m 
%
% uicontrol callback function used to let the user either repostion axial
% plane or select MRSI spec from the popup menu in the spm graphics. This method sets other callback methods to
% allow theese functions to occur. For more info about callback functions
% are found here 
% https://www.mathworks.com/help/matlab/creating_plots/callback-definition.html
%
% USAGE:
% Can't use (callback function)
%
% INPUT:
% src = returned updated uicontrol element
%
function selection(src, ~)
    %get the spm global variable
    global st
    %get the axis object for the axial plane
    axis = st.vols{1}.ax{1}.ax;
    
    if(src.Value == 2)
       %If the popup menu is set to select_spec
       %set the callback function of click on axis to the method
       %select_spec.m
       set(axis, 'ButtonDownFcn', @select_spec);
       %position new axis here with properties
       ax = axes('Visible','off', 'Parent',st.fig, ...
        'YDir','normal', 'Units','Pixels', 'Box','on',...
        'Position',[340 445 250 220],...
        'Units','normalized',...
        'Visible','on', 'Tag', 'csi_axis');
        
    
    else
        %If the popup menu is set to reposition
        %set click on axial plane to reposition
        set(axis, 'ButtonDownFcn', @repos_start);
        %find the axis used to plot MRSI specs
        ax = findobj(gcf, 'Tag', 'csi_axis');
        %delete the specs
        delete(ax);
    end
    
    %callback method used for repositioning
    function repos_start(~, ~)
        
        %checking if the click is a right click
        if ~strcmpi(get(gcbf,'SelectionType'),'alt')
            %set callback functions for movement of mouse and finished
            %clicking. spm_orthviews('reposition') is for movement of mouse
            % and repos_end for end of click.
            set(gcbf,'windowbuttonmotionfcn',@()(spm_orthviews('reposition')), 'windowbuttonupfcn',@repos_end);
            %call reposition in spm_orthviews
            spm_orthviews('reposition');
        end
    end
   
    function repos_end(~, ~)
        %remove callback functions for movement of mouse and wiindows
        %button up.
        set(gcbf,'windowbuttonmotionfcn','', 'windowbuttonupfcn','');
    end       
end