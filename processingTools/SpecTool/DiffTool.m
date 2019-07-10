function varargout = DiffTool(varargin)
%DiffTool.m
%Jamie Near, McGill University 2014.  
%
% USAGE:
% DiffTool(in1,in2,ppmmin,ppmmax);
% 
% DESCRIPTION:
% Graphical User Interface for Subtraction of two spectra, with manual 
% adjustment of the relative phase and frequency of the two spectra.  
% Difference spectrum is shown in real time.  
% 
% INPUTS:
% in1        = first input dataset in matlab structure format.
% in2        = second input dataset in matlab structure format.
% ppmmin     = minimum of ppm frequency range.
% ppmmax     = maximum of ppm frequency range.
%
%DIFFTOOL M-file for DiffTool.fig
%      DIFFTOOL, by itself, creates a new DIFFTOOL or raises the existing
%      singleton*.
%
%      H = DIFFTOOL returns the handle to a new DIFFTOOL or the handle to
%      the existing singleton*.
%
%      DIFFTOOL('Property','Value',...) creates a new DIFFTOOL using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to DiffTool_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      DIFFTOOL('CALLBACK') and DIFFTOOL('CALLBACK',hObject,...) call the
%      local function named CALLBACK in DIFFTOOL.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DiffTool

% Last Modified by GUIDE v2.5 16-Jul-2014 12:48:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DiffTool_OpeningFcn, ...
                   'gui_OutputFcn',  @DiffTool_OutputFcn, ...
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


% --- Executes just before DiffTool is made visible.
function DiffTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

%JAMIE NEAR ADDED SPECTROSCOPY DATA CODE

handles.in1=varargin{1};
handles.in2=varargin{2};
handles.ppmmin=varargin{3};
handles.ppmmax=varargin{4};

%Magnetic field 3 Tesla
Bo=handles.in1.Bo;

n=length(handles.in1.t);

handles.orig_Specs1=handles.in1.specs;
handles.orig_Specs2=handles.in2.specs;
handles.current_Specs2=handles.in2.specs;
handles.diff=op_subtractScans(handles.in1,handles.in2);
handles.orig_diffSpecs=handles.diff.specs;
handles.current_diffSpecs=handles.diff.specs;

%obtains the zero order phase value from the phas0Slider component
phas0Value = 0;
phas1Value = 0;
freqValue = 0;
ampValue = 100;
DCValue = 0;
 
handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;



labsiz=14;
handles.labsiz=labsiz;
tiksiz=12;
handles.tiksiz=tiksiz;
%linwidth=1.5;
%mksiz=12;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);


%END OF SECTION ADDED BY JAMIE NEAR

% Choose default command line output for DiffTool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DiffTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = DiffTool_OutputFcn(hObject, eventdata, handles)
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
function phas0Slider_Callback(hObject, eventdata, handles)
% hObject    handle to phas0Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%obtains the zero order phase value from the phas0Slider component
phas0Value = get(handles.phas0Slider,'Value');
phas1Value = get(handles.phas1Slider,'Value');
freqValue = get(handles.freqSlider,'Value');
ampValue = get(handles.ampSlider,'Value');
DCValue = get(handles.DCSlider,'Value');
 
%puts the slider value into the edit text component
set(handles.phas0Slider_editText,'String', num2str(phas0Value));

handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

% Update handles structure
guidata(hObject, handles);






% --- Executes during object creation, after setting all properties.
function phas0Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phas0Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end






% --- Executes on key press with focus on phas0Slider and none of its controls.
function phas0Slider_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to phas0Slider (see GCBO)
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



function phas0Slider_editText_Callback(hObject, eventdata, handles)
% hObject    handle to phas0Slider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phas0Slider_editText as text
%        str2double(get(hObject,'String')) returns contents of phas0Slider_editText as a double

%get the string for the editText component
phas0Value = get(handles.phas0Slider_editText,'String');
phas1Value = get(handles.phas1Slider,'String');
freqValue = get(handles.freqSlider,'Value');
ampValue = get(handles.ampSlider,'Value');
DCValue = get(handles.DCSlider,'Value');
 
%convert from string to number if possible, otherwise returns empty
phas0Value = str2num(phas0Value);

 
%if user inputs something is not a number, or if the input is less than 0
%then the slider value defaults to 0
if (isempty(phas0Value) || phas0Value < -180 || phas0Value > 180)
    set(handles.phas0Slider,'Value',0);
    set(handles.phas0Slider_editText,'String','0');
else
    set(handles.phas0Slider,'Value',phas0Value);
end
handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

% Update handles structure
guidata(hObject, handles);







% --- Executes during object creation, after setting all properties.
function phas0Slider_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phas0Slider_editText (see GCBO)
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
phas0Value = get(handles.phas0Slider,'Value');
phas1Value = get(handles.phas1Slider,'Value');
freqValue = get(handles.freqSlider,'Value');
ampValue = get(handles.ampSlider,'Value');
DCValue = get(handles.DCSlider,'Value');
 
%puts the slider value into the edit text component
set(handles.freqSlider_editText,'String', num2str(freqValue));

% Update handles structure
guidata(hObject, handles);


handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

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
phas0Value = get(handles.phas0Slider,'Value');
phas1Value = get(handles.phas1Slider,'Value');
freqValue = get(handles.freqSlider_editText,'String');
ampValue = get(handles.ampSlider,'Value');
DCValue = get(handles.DCSlider,'Value');
 
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

handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

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


% --- Executes on key press with focus on phas0Slider_editText and none of its controls.
function phas0Slider_editText_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to phas0Slider_editText (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function phas1Slider_Callback(hObject, eventdata, handles)
% hObject    handle to phas1Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%obtains the zero order phase value from the phas0Slider component
phas0Value = get(handles.phas0Slider,'Value');
phas1Value = get(handles.phas1Slider,'Value');
freqValue = get(handles.freqSlider,'Value');
ampValue = get(handles.ampSlider,'Value');
DCValue = get(handles.DCSlider,'Value');
 
%puts the slider value into the edit text component
set(handles.phas1Slider_editText,'String', num2str(phas1Value));

handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function phas1Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phas1Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function phas1Slider_editText_Callback(hObject, eventdata, handles)
% hObject    handle to phas1Slider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phas1Slider_editText as text
%        str2double(get(hObject,'String')) returns contents of phas1Slider_editText as a double

%get the string for the editText component
phas0Value = get(handles.phas0Slider,'String');
phas1Value = get(handles.phas1Slider_editText,'String');
freqValue = get(handles.freqSlider,'Value');
ampValue = get(handles.ampSlider,'Value');
DCValue = get(handles.DCSlider,'Value');
 
%convert from string to number if possible, otherwise returns empty
phas1Value = str2num(phas1Value);

 
%if user inputs something is not a number, or if the input is less than 0
%then the slider value defaults to 0
if (isempty(phas1Value) || phas1Value < -1 || phas1Value > 1)
    set(handles.phas1Slider,'Value',0);
    set(handles.phas1Slider_editText,'String','0');
else
    set(handles.phas1Slider,'Value',phas1Value);
end
handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

% Update handles structure
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function phas1Slider_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phas1Slider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function ampSlider_Callback(hObject, eventdata, handles)
% hObject    handle to ampSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%obtains the zero order phase value from the phas0Slider component
phas0Value = get(handles.phas0Slider,'Value');
phas1Value = get(handles.phas1Slider,'Value');
freqValue = get(handles.freqSlider,'Value');
ampValue = get(handles.ampSlider,'Value');
DCValue = get(handles.DCSlider,'Value');
 
%puts the slider value into the edit text component
set(handles.ampSlider_editText,'String', num2str(ampValue));

handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ampSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ampSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function ampSlider_editText_Callback(hObject, eventdata, handles)
% hObject    handle to ampSlider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ampSlider_editText as text
%        str2double(get(hObject,'String')) returns contents of ampSlider_editText as a double

%get the string for the editText component
phas0Value = get(handles.phas0Slider,'String');
phas1Value = get(handles.phas1Slider,'String');
freqValue = get(handles.freqSlider,'Value');
ampValue = get(handles.ampSlider_editText,'Value');
DCValue = get(handles.DCSlider,'Value');
 
%convert from string to number if possible, otherwise returns empty
ampValue = str2num(ampValue);

 
%if user inputs something is not a number, or if the input is less than 0
%then the slider value defaults to 0
if (isempty(ampValue) || ampValue < 0 || ampValue > 200)
    set(handles.ampSlider,'Value',100);
    set(handles.ampSlider_editText,'String','100');
else
    set(handles.ampSlider,'Value',ampValue);
end
handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function ampSlider_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ampSlider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function DCSlider_Callback(hObject, eventdata, handles)
% hObject    handle to DCSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

%obtains the zero order phase value from the phas0Slider component
phas0Value = get(handles.phas0Slider,'Value');
phas1Value = get(handles.phas1Slider,'Value');
freqValue = get(handles.freqSlider,'Value');
ampValue = get(handles.ampSlider,'Value');
DCValue = get(handles.DCSlider,'Value');
 
%puts the slider value into the edit text component
set(handles.DCSlider_editText,'String', num2str(DCValue));

handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function DCSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DCSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function DCSlider_editText_Callback(hObject, eventdata, handles)
% hObject    handle to DCSlider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DCSlider_editText as text
%        str2double(get(hObject,'String')) returns contents of DCSlider_editText as a double

%get the string for the editText component
phas0Value = get(handles.phas0Slider,'String');
phas1Value = get(handles.phas1Slider,'String');
freqValue = get(handles.freqSlider,'Value');
ampValue = get(handles.ampSlider,'Value');
DCValue = get(handles.DCSlider_editText,'Value');
 
%convert from string to number if possible, otherwise returns empty
DCValue = str2num(DCValue);

 
%if user inputs something is not a number, or if the input is less than 0
%then the slider value defaults to 0
if (isempty(DCValue) || DCValue < -500 || DCValue > 500)
    set(handles.DCSlider,'Value',0);
    set(handles.DCSlider_editText,'String','0');
else
    set(handles.DCSlider,'Value',DCValue);
end
handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;
handles.freqValue=freqValue;
handles.ampValue=ampValue;
handles.DCValue=DCValue;

%Create current_timeDomainData and current freqDomainData
in2_phas=op_addphase(handles.in2,handles.phas0Value,handles.phas1Value);
in2_phas_freq=op_freqshift(in2_phas,handles.freqValue);
in2_phas_freq_amp=op_ampScale(in2_phas_freq,handles.ampValue/100);
%in2_phas_freq_amp_dc=op_dccorr(in2_phas_freq_amp,'v',handles.DCValue*max(real(handles.in1.specs))/5000);
in2_phas_freq_amp_dc=op_filter2(in2_phas_freq_amp,'g',handles.DCValue);
handles.current_Specs2=in2_phas_freq_amp_dc.specs;
handles.diff=op_subtractScans(handles.in1,in2_phas_freq_amp_dc);
handles.current_diffSpecs=handles.diff.specs;

% %MAKE SUBSPECS PLOTS
axes(handles.axes2);
plot(handles.in1.ppm,real(handles.orig_Specs1),'b-',handles.in2.ppm,real(handles.current_Specs2),'r-',handles.diff.ppm,real(handles.current_diffSpecs)*5,'g-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function DCSlider_editText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DCSlider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
