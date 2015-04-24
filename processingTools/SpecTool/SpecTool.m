function varargout = SpecTool(varargin)
%SpecTool.m
%Jamie Near, McGill University 2014.  
%
% USAGE:
% SpecTool(in,tmax,ppmmin,ppmmax);
% 
% DESCRIPTION:
% Graphical User Interface for manual adjustment of zero- and first-order
% phase of an MR spectrum.  
% 
% INPUTS:
% in         = input dataset in matlab structure format.
% tmax       = maximum timepoint (s) for fid viewing
% ppmmin     = minimum of ppm frequency range.
% ppmmax     = maximum of ppm frequency range.
%
% SPECTOOL M-file for SpecTool.fig
%      SPECTOOL, by itself, creates a new SPECTOOL or raises the existing
%      singleton*.
%
%      H = SPECTOOL returns the handle to a new SPECTOOL or the handle to
%      the existing singleton*.
%
%      SPECTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECTOOL.M with the given input arguments.
%
%      SPECTOOL('Property','Value',...) creates a new SPECTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SpecTool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SpecTool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SpecTool

% Last Modified by GUIDE v2.5 26-Nov-2010 13:30:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SpecTool_OpeningFcn, ...
                   'gui_OutputFcn',  @SpecTool_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before SpecTool is made visible.
function SpecTool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SpecTool (see VARARGIN)


%JAMIE NEAR ADDED SYNTHETIC SPECTROSCOPY DATA CODE

handles.in=varargin{1};
handles.tmax=varargin{2};
if length(varargin)>3
    handles.ppmmin=varargin{3};
    handles.ppmmax=varargin{4};
else
    handles.ppmmin=min(handles.in.ppm);
    handles.ppmmax=max(handles.in.ppm);
end

handles.ppmmin
handles.ppmmax

%Magnetic field 3 Tesla
Bo=handles.in.Bo;

n=length(handles.in.t);

handles.orig_timeDomainData=handles.in.fids;
handles.current_timeDomainData=handles.in.fids;
handles.orig_freqDomainData=handles.in.specs;
handles.current_freqDomainData=handles.in.specs;

%obtains the zero order phase value from the phas0Slider component
phas0Value = 0;
phas1Value = 0;
 
handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;



labsiz=14;
handles.labsiz=labsiz;
tiksiz=12;
handles.tiksiz=tiksiz;
%linwidth=1.5;
%mksiz=12;

% %MAKE TIME DOMAIN PLOTS
axes(handles.axes2);
plot(handles.in.t,real(handles.current_timeDomainData),'b-');
xlim([0 handles.tmax]);
%set(gca,'FontSize',handles.tiksiz)
%set(gcf,'Color','w');
%box off;
%xlabel('Time (s)','FontSize',handles.labsiz);
%ylabel('Signal Intensity','FontSize',handles.labsiz);


%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_freqDomainData),'b-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'XDir','reverse');
%set(gca,'FontSize',tiksiz);
%set(gcf,'Color','w');
%box off;
%xlabel('Frequency (ppm)','FontSize',labsiz);
%ylabel('Spectral Intensity','FontSize',labsiz);


%END OF SECTION ADDED BY JAMIE NEAR

% Choose default command line output for SpecTool
handles.output=hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SpecTool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SpecTool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





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
 
%puts the slider value into the edit text component
set(handles.phas0Slider_editText,'String', num2str(phas0Value));

handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;

%Create current_timeDomainData and current freqDomainData
in_ph=op_addphase(handles.in,handles.phas0Value,handles.phas1Value);
handles.current_timeDomainData=in_ph.fids;
handles.current_freqDomainData=in_ph.specs;

% %MAKE TIME DOMAIN PLOTS
axes(handles.axes2);
plot(handles.in.t,real(handles.current_timeDomainData),'b-');
xlim([0 handles.tmax]);
% set(gca,'FontSize',handles.tiksiz)
% set(gcf,'Color','w');
% box off;
% xlabel('Time (s)','FontSize',handles.labsiz);
% ylabel('Signal Intensity','FontSize',handles.labsiz);

%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_freqDomainData),'b-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'XDir','reverse');
% set(gca,'FontSize',handles.tiksiz);
% set(gcf,'Color','w');
% box off;
% xlabel('Frequency (ppm)','FontSize',handles.labsiz);
% ylabel('Spectral Intensity','FontSize',handles.labsiz);

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






function phas0Slider_editText_Callback(hObject, eventdata, handles)
% hObject    handle to phas0Slider_editText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phas0Slider_editText as text
%        str2double(get(hObject,'String')) returns contents of phas0Slider_editText as a double

%get the string for the editText component
phas0Value = get(handles.phas0Slider_editText,'String');
phas1Value = get(handles.phas1Slider,'Value');
 
%convert from string to number if possible, otherwise returns empty
phas0Value = str2num(phas0Value);

 
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(phas0Value) || phas0Value < -270 || phas0Value > 270)
    set(handles.phas0Slider,'Value',0);
    set(handles.phas0Slider_editText,'String','0');
else
    set(handles.phas0Slider,'Value',phas0Value);
end
handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;

%Create current_timeDomainData and current freqDomainData
in_ph=op_addphase(handles.in,handles.phas0Value,handles.phas1Value);
handles.current_timeDomainData=in_ph.fids;
handles.current_freqDomainData=in_ph.specs;

% %MAKE TIME DOMAIN PLOTS
axes(handles.axes2);
plot(handles.in.t,real(handles.current_timeDomainData),'b-');
xlim([0 handles.tmax]);
% set(gca,'FontSize',handles.tiksiz)
% set(gcf,'Color','w');
% box off;
% xlabel('Time (s)','FontSize',handles.labsiz);
% ylabel('Signal Intensity','FontSize',handles.labsiz);

%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_freqDomainData),'b-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'XDir','reverse');
% set(gca,'FontSize',handles.tiksiz);
% set(gcf,'Color','w');
% box off;
% xlabel('Frequency (ppm)','FontSize',handles.labsiz);
% ylabel('Spectral Intensity','FontSize',handles.labsiz);

% Update handles structure
guidata(hObject,handles);







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
function phas1Slider_Callback(hObject, eventdata, handles)
% hObject    handle to phas1Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider



%obtains the slider value from the slider component
phas0Value = get(handles.phas0Slider,'Value');
phas1Value = get(handles.phas1Slider,'Value');
 
%puts the slider value into the edit text component
set(handles.phas1Slider_editText,'String', num2str(phas1Value));

% Update handles structure
guidata(hObject, handles);


handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;

%Create current_timeDomainData and current freqDomainData
in_ph=op_addphase(handles.in,handles.phas0Value,handles.phas1Value);
handles.current_timeDomainData=in_ph.fids;
handles.current_freqDomainData=in_ph.specs;

% %MAKE TIME DOMAIN PLOTS
axes(handles.axes2);
plot(handles.in.t,real(handles.current_timeDomainData),'b-');
xlim([0 handles.tmax]);
% set(gca,'FontSize',handles.tiksiz)
% set(gcf,'Color','w');
% box off;
% xlabel('Time (s)','FontSize',handles.labsiz);
% ylabel('Signal Intensity','FontSize',handles.labsiz);

%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_freqDomainData),'b-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'XDir','reverse');
% set(gca,'FontSize',handles.tiksiz);
% set(gcf,'Color','w');
% box off;
% xlabel('Frequency (ppm)','FontSize',handles.labsiz);
% ylabel('Spectral Intensity','FontSize',handles.labsiz);

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
%        str2double(get(hObject,'String')) returns contents of
%        phas1Slider_editText as a double

%get the string for the editText component
phas0Value = get(handles.phas0Slider,'Value');
phas1Value = get(handles.phas1Slider_editText,'String');
 
%convert from string to number if possible, otherwise returns empty
phas1Value = str2num(phas1Value);
 
%if user inputs something is not a number, or if the input is less than 0
%or greater than 100, then the slider value defaults to 0
if (isempty(phas1Value) || phas1Value < -0.05 || phas1Value > 0.05   )
    phas1Value=0;
    set(handles.phas1Slider,'Value',0);
    set(handles.phas1Slider_editText,'String','0');
else
    set(handles.phas1Slider,'Value',phas1Value);
end

handles.phas0Value=phas0Value;
handles.phas1Value=phas1Value;

%Create current_timeDomainData and current freqDomainData
in_ph=op_addphase(handles.in,handles.phas0Value,handles.phas1Value);
handles.current_timeDomainData=in_ph.fids;
handles.current_freqDomainData=in_ph.specs;

% %MAKE TIME DOMAIN PLOTS
axes(handles.axes2);
plot(handles.in.t,real(handles.current_timeDomainData),'b-');
xlim([0 handles.tmax]);
% set(gca,'FontSize',handles.tiksiz)
% set(gcf,'Color','w');
% box off;
% xlabel('Time (s)','FontSize',handles.labsiz);
% ylabel('Signal Intensity','FontSize',handles.labsiz);

%MAKE FREQENCY DOMAIN PLOTS
axes(handles.axes3);
plot(handles.in.ppm,real(handles.current_freqDomainData),'b-');
xlim([handles.ppmmin handles.ppmmax]);
%set(gca,'XDir','reverse');
% set(gca,'FontSize',handles.tiksiz);
% set(gcf,'Color','w');
% box off;
% xlabel('Frequency (ppm)','FontSize',handles.labsiz);
% ylabel('Spectral Intensity','FontSize',handles.labsiz);


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




function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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


% --- Executes on key press with focus on phas0Slider and none of its controls.
function phas0Slider_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to phas0Slider (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
