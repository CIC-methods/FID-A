function varargout = subSpecTool(varargin)
%subSpecTool.m
%Jamie Near, McGill University 2014.
%
% USAGE:
% subSpecTool(in,ppmmin,ppmmax);
% 
% DESCRIPTION:
% Graphical User Interface for manual adjustment of relative frequency and
% phase of edit-on and edit-off spectra for a MEGA-PRESS dataset.
% 
% INPUTS:
% in         = input dataset in matlab structure format.
% ppmmin     = minimum of ppm frequency range.
% ppmmax     = maximum of ppm frequency range.

%SUBSPECTOOL M-file for subSpecTool.fig
%      SUBSPECTOOL, by itself, creates a new SUBSPECTOOL or raises the existing
%      singleton*.
%
%      H = SUBSPECTOOL returns the handle to a new SUBSPECTOOL or the handle to
%      the existing singleton*.
%
%      SUBSPECTOOL('Property','Value',...) creates a new SUBSPECTOOL using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to subSpecTool_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SUBSPECTOOL('CALLBACK') and SUBSPECTOOL('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SUBSPECTOOL.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help subSpecTool

% Last Modified by GUIDE v2.5 03-Apr-2013 15:35:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @subSpecTool_OpeningFcn, ...
                   'gui_OutputFcn',  @subSpecTool_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
   gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before subSpecTool is made visible.
function subSpecTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

%JAMIE NEAR ADDED SPECTROSCOPY DATA CODE

handles.in=varargin{1};
handles.ppmmin=varargin{2};
handles.ppmmax=varargin{3};

%Magnetic field 3 Tesla
Bo=handles.in.Bo;

n=length(handles.in.t);

handles.orig_subSpecs=handles.in.specs;
handles.current_subSpecs=handles.in.specs;
handles.diff=op_combinesubspecs(handles.in,'diff');
handles.orig_diffSpecs=handles.diff.specs;
handles.current_diffSpecs=handles.diff.specs;

%obtains the zero order phase value from the phas0Slider component
phasValue = 0;
freqValue = 0;
 
handles.phasValue=phasValue;
handles.freqValue=freqValue;



labsiz=14;
handles.labsiz=labsiz;
tiksiz=12;
handles.tiksiz=tiksiz;
%linwidth=1.5;
%mksiz=12;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in.ppm,real(handles.current_subSpecs(:,1)),'b-',handles.in.ppm,real(addphase(handles.current_subSpecs(:,2),180)),'r-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);


%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_diffSpecs),'b-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',tiksiz);
%set(gcf,'Color','w');
%box off;
%xlabel('Frequency (ppm)','FontSize',labsiz);
%ylabel('Spectral Intensity','FontSize',labsiz);


%END OF SECTION ADDED BY JAMIE NEAR

% Choose default command line output for subSpecTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes subSpecTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = subSpecTool_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes on slider movement.
function phasSlider_Callback(hObject, eventdata, handles)
% hObject    handle to phasSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%obtains the zero order phase value from the phas0Slider component
phasValue = get(handles.phasSlider,'Value');
freqValue = get(handles.freqSlider,'Value');
 
%puts the slider value into the edit text component
set(handles.phasSlider_editText,'String', num2str(phasValue));

handles.phasValue=phasValue;
handles.freqValue=freqValue;

%Create current_timeDomainData and current freqDomainData
in_phas=op_addphaseSubspec(handles.in,handles.phasValue);
in_phas_freq=op_freqshiftSubspec(in_phas,handles.freqValue);
handles.current_subSpecs=in_phas_freq.specs;
handles.diff=op_combinesubspecs(in_phas_freq,'diff');
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in.ppm,real(handles.current_subSpecs(:,1)),'b-',handles.in.ppm,real(addphase(handles.current_subSpecs(:,2),180)),'r-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);


%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_diffSpecs),'b-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',tiksiz);
%set(gcf,'Color','w');
%box off;
%xlabel('Frequency (ppm)','FontSize',labsiz);
%ylabel('Spectral Intensity','FontSize',labsiz);

% Update handles structure
guidata(hObject, handles);






% --- Executes during object creation, after setting all properties.
function phasSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phasSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end






% --- Executes on key press with focus on phasSlider and none of its controls.
function phasSlider_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to phasSlider (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function phasSlider_editText_Callback(hObject, eventdata, handles)
% hObject    handle to phasSlider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phasSlider_editText as text
%        str2double(get(hObject,'String')) returns contents of phasSlider_editText as a double

%get the string for the editText component
phasValue = get(handles.phasSlider_editText,'String');
freqValue = get(handles.freqSlider,'Value');
 
%convert from string to number if possible, otherwise returns empty
phasValue = str2num(phasValue);

 
%if user inputs something is not a number, or if the input is less than 0
%then the slider value defaults to 0
if (isempty(phasValue) || phasValue < -180 || phasValue > 180)
    set(handles.phasSlider,'Value',0);
    set(handles.phasSlider_editText,'String','0');
else
    set(handles.phasSlider,'Value',phasValue);
end
handles.phasValue=phasValue;
handles.freqValue=freqValue;

%Create current_timeDomainData and current freqDomainData
in_phas=op_addphaseSubspec(handles.in,handles.phasValue);
in_phas_freq=op_freqshiftSubspec(in_phas,handles.freqValue);
handles.current_subSpecs=in_phas_freq.specs;
handles.diff=op_combinesubspecs(in_phas_freq,'diff');
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in.ppm,real(handles.current_subSpecs(:,1)),'b-',handles.in.ppm,real(addphase(handles.current_subSpecs(:,2),180)),'r-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);


%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_diffSpecs),'b-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',tiksiz);
%set(gcf,'Color','w');
%box off;
%xlabel('Frequency (ppm)','FontSize',labsiz);
%ylabel('Spectral Intensity','FontSize',labsiz);

% Update handles structure
guidata(hObject, handles);







% --- Executes during object creation, after setting all properties.
function phasSlider_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phasSlider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on slider movement.
function freqSlider_Callback(hObject, eventdata, handles)
% hObject    handle to freqSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


%obtains the slider value from the slider component
phasValue = get(handles.phasSlider,'Value');
freqValue = get(handles.freqSlider,'Value');
 
%puts the slider value into the edit text component
set(handles.freqSlider_editText,'String', num2str(freqValue));

% Update handles structure
guidata(hObject, handles);


handles.phasValue=phasValue;
handles.freqValue=freqValue;

%Create current_timeDomainData and current freqDomainData
in_phas=op_addphaseSubspec(handles.in,handles.phasValue);
in_phas_freq=op_freqshiftSubspec(in_phas,handles.freqValue);
handles.current_subSpecs=in_phas_freq.specs;
handles.diff=op_combinesubspecs(in_phas_freq,'diff');
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in.ppm,real(handles.current_subSpecs(:,1)),'b-',handles.in.ppm,real(addphase(handles.current_subSpecs(:,2),180)),'r-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);


%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_diffSpecs),'b-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',tiksiz);
%set(gcf,'Color','w');
%box off;
%xlabel('Frequency (ppm)','FontSize',labsiz);
%ylabel('Spectral Intensity','FontSize',labsiz);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function freqSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function freqSlider_editText_Callback(hObject, eventdata, handles)
% hObject    handle to freqSlider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of freqSlider_editText as text
%        str2double(get(hObject,'String')) returns contents of freqSlider_editText as a double


%get the string for the editText component
phasValue = get(handles.phasSlider,'Value');
freqValue = get(handles.freqSlider_editText,'String');
 
%convert from string to number if possible, otherwise returns empty
freqValue = str2num(freqValue);
 
%if user inputs something is not a number, or if the input is less than
%-100 or greater than 100, then the slider value defaults to 0
if (isempty(freqValue) || freqValue < -100 || freqValue > 100   )
    freqValue=0;
    set(handles.freqSlider,'Value',0);
    set(handles.freqSlider_editText,'String','0');
else
    set(handles.freqSlider,'Value',freqValue);
end

handles.phasValue=phasValue;
handles.freqValue=freqValue;

%Create current_timeDomainData and current freqDomainData
in_phas=op_addphaseSubspec(handles.in,handles.phasValue);
in_phas_freq=op_freqshiftSubspec(in_phas,handles.freqValue);
handles.current_subSpecs=in_phas_freq.specs;
handles.diff=op_combinesubspecs(in_phas_freq,'diff');
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in.ppm,real(handles.current_subSpecs(:,1)),'b-',handles.in.ppm,real(addphase(handles.current_subSpecs(:,2),180)),'r-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',handles.tiksiz);
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);


%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_diffSpecs),'b-','LineWidth',1.2);
xlim([handles.ppmmin handles.ppmmax]);
set(gca,'XDir','reverse');
%set(gca,'FontSize',tiksiz);
%set(gcf,'Color','w');
%box off;
%xlabel('Frequency (ppm)','FontSize',labsiz);
%ylabel('Spectral Intensity','FontSize',labsiz);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function freqSlider_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to freqSlider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on phasSlider_editText and none of its controls.
function phasSlider_editText_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to phasSlider_editText (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
