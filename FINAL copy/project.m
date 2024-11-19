
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Auhter: Sarina Ziraksima %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout = project(varargin)
% PROJECT MATLAB code for project.fig
%      PROJECT, by itself, creates a new PROJECT or raises the existing
%      singleton*.
%
%      H = PROJECT returns the handle to a new PROJECT or the handle to
%      the existing singleton*.
%
%      PROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECT.M with the given input arguments.
%
%      PROJECT('Property','Value',...) creates a new PROJECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before project_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to project_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help project

% Last Modified by GUIDE v2.5 17-Jun-2024 02:44:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @project_OpeningFcn, ...
                   'gui_OutputFcn',  @project_OutputFcn, ...
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


% --- Executes just before project is made visible.
function project_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

            %set(handles.Apass_edit,'Enable','off');
            %set(handles.Astop_edit,'Enable','off');
            %set(handles.Fpass_edit,'Enable','off');
            %set(handles.Fstop_edit,'Enable','off');
            
            %set(handles.fn_edit,'Enable','off');
            %set(handles.N_edit,'Enable','off');
            %output_img=imread("blank.png");
            
            %imshow(output_img)
            
            set(handles.table1,'Visible','off');
            set(handles.text21,'Visible','off');
            set(handles.text24,'Visible','off');
            set(handles.text28,'Visible','off');
            set(handles.stages_pop,'Visible','off');
            
            

% varargin   command line arguments to project (see VARARGIN)

% Choose default command line output for project
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes project wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = project_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



Apass = str2num(get(handles.Apass_edit,'string'));
Astop = str2num(get(handles.Astop_edit,'string'));
handles.N = str2num(get(handles.N_edit,'string'));
fn = str2num(get(handles.fn_edit,'string'));
Fpass = str2num(get(handles.Fpass_edit,'string'));
Fstop = str2num(get(handles.Fstop_edit,'string'));
n = str2num(get(handles.n_edit,'string'));
gain = str2num(get(handles.gain_edit,'string'));



switch (get(handles.opamps_pop,'Value'))
    case 1  % costumize      [Ri Ro wu gain cost]
        opamps = [1e6 400 6e5 10e6 100e3 0];
    case 2  % LM324
        opamps = [10e6 300 1.2e6 100e3 32800];
    case 3  % UA741
        opamps = [2e6 75 1e6 200e3 43700];
    case 4  % TL072 JFET
        opamps = [100e6 125 1e6 1.7783e+06 38000];  
    case 5  % TL322C
        opamps = [1e6 75 1e6 200e3 132000];
    case 6  % LM2902N
        opamps = [10e6 300 1.2e6 100e3 50000];
    
        
end
% 'low' | 'bandpass' | 'high' | 'stop'
switch (get(handles.ftype_pop,'Value'))
    case 1        
        handles.ftype = 'low';
    case 2
        handles.ftype = 'bandpass';
    case 3
        handles.ftype = 'high';
    case 4
        handles.ftype = 'stop';  
       
end

% 'butter' | 'cheby1' | 'cheby2' | 'ellip'
switch (get(handles.methoud_solution_pop,'Value'))
    case 1        
        handles.methodSolution = 'butter';
    case 2
        handles.methodSolution = 'cheby1';
    case 3
        handles.methodSolution = 'cheby2';
    case 4
        handles.methodSolution = 'ellip';  
       
end
%'sallen key' | 'MFB' | 'akerberg' | 'thomas 1' | 'thomas 2' | 'optimize'
switch (get(handles.topology_type_pop,'Value'))
    case 1        
        handles.topology_type = 'sallen key';
    case 2
        handles.topology_type = 'MFB';
    case 3
        handles.topology_type = 'akerberg';
    case 4
        handles.topology_type = 'thomas 1';  
    case 5
        handles.topology_type = 'thomas 2';
    case 6
        handles.topology_type = 'optimize';
        
end

% resister and capaciture values
switch (get(handles.r_c_val_pop,'Value'))
    case 1 % costumize        
        r_c_val = 0;
    case 2 % E24
        r_c_val = 1;
    %case 3
    %    topology_type = 'akerberg';
    %case 4
    %    topology_type = 'thomas 1';  
    %case 5
    %    topology_type = 'thomas 2';
    %case 6
    %    topology_type = 'optimize';
        
end


%mode




FitnessLimit_Data = str2num(get(handles.FitnessLimit_Data_edit,'string'));
MaxGenerations_Data = str2num(get(handles.MaxGenerations_Data_edit,'string'));



[handles.errorr,handles.final_results]=main_pj(gain,r_c_val,fn,Fpass,Apass,Fstop,Astop,handles.methodSolution,handles.ftype,handles.N,handles.topology_type,n,FitnessLimit_Data,MaxGenerations_Data,opamps)
set(handles.table1,'Visible','on');
set(handles.text21,'Visible','on');
set(handles.text24,'Visible','on');
set(handles.text28,'Visible','on');
set(handles.stages_pop,'Visible','on');
set(handles.table1,'data',handles.final_results(1,:));
set(handles.stages_pop,'String',{[1:length(handles.final_results(:,1))]});
set(handles.errorr_text,'Visible','on');
set(handles.errorr_text,'String',handles.errorr);





%handles.output = hObject;
guidata(hObject,handles);

function Apass_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Apass_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Apass_edit as text
%        str2double(get(hObject,'String')) returns contents of Apass_edit as a double


% --- Executes during object creation, after setting all properties.
function Apass_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Apass_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Astop_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Astop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Astop_edit as text
%        str2double(get(hObject,'String')) returns contents of Astop_edit as a double


% --- Executes during object creation, after setting all properties.
function Astop_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Astop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function N_edit_Callback(hObject, eventdata, handles)
% hObject    handle to N_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_edit as text
%        str2double(get(hObject,'String')) returns contents of N_edit as a double


% --- Executes during object creation, after setting all properties.
function N_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function fn_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fn_edit as text
%        str2double(get(hObject,'String')) returns contents of fn_edit as a double


% --- Executes during object creation, after setting all properties.
function fn_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fn_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Fpass_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Fpass_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fpass_edit as text
%        str2double(get(hObject,'String')) returns contents of Fpass_edit as a double


% --- Executes during object creation, after setting all properties.
function Fpass_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fpass_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Fstop_edit_Callback(hObject, eventdata, handles)
% hObject    handle to Fstop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fstop_edit as text
%        str2double(get(hObject,'String')) returns contents of Fstop_edit as a double


% --- Executes during object creation, after setting all properties.
function Fstop_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fstop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function n_edit_Callback(hObject, eventdata, handles)
% hObject    handle to n_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of n_edit as text
%        str2double(get(hObject,'String')) returns contents of n_edit as a double


% --- Executes during object creation, after setting all properties.
function n_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to n_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function FitnessLimit_Data_edit_Callback(hObject, eventdata, handles)
% hObject    handle to FitnessLimit_Data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of FitnessLimit_Data_edit as text
%        str2double(get(hObject,'String')) returns contents of FitnessLimit_Data_edit as a double


% --- Executes during object creation, after setting all properties.
function FitnessLimit_Data_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FitnessLimit_Data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxGenerations_Data_edit_Callback(hObject, eventdata, handles)
% hObject    handle to MaxGenerations_Data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxGenerations_Data_edit as text
%        str2double(get(hObject,'String')) returns contents of MaxGenerations_Data_edit as a double


% --- Executes during object creation, after setting all properties.
function MaxGenerations_Data_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxGenerations_Data_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end













% --- Executes on selection change in topology_type_pop.
function topology_type_pop_Callback(hObject, eventdata, handles)
% hObject    handle to topology_type_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns topology_type_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from topology_type_pop


% --- Executes during object creation, after setting all properties.
function topology_type_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to topology_type_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in opamps_pop.
function opamps_pop_Callback(hObject, eventdata, handles)
% hObject    handle to opamps_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns opamps_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from opamps_pop


% --- Executes during object creation, after setting all properties.
function opamps_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to opamps_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ftype_pop.
function ftype_pop_Callback(hObject, eventdata, handles)
% hObject    handle to ftype_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ftype_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ftype_pop


% --- Executes during object creation, after setting all properties.
function ftype_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ftype_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in stages_pop.
function stages_pop_Callback(hObject, eventdata, handles)
% hObject    handle to stages_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

     
 set(handles.table1,'data',handles.final_results(get(handles.stages_pop,'Value'),:));
 tiledlayout(2,1)
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%get(handles.stages_pop,'Value'),:)





if(handles.ftype(1)=='l')
    handles.x2='L';
end

if(handles.ftype(1)=='h') 
    handles.x2='H';
end

if(handles.ftype(1)=='b')
    handles.x2='B';
end


if(handles.topology_type(1)=='s')
    handles.x1='S_';
end

if(handles.topology_type(1)=='M') 
    handles.x1='M_';
end

if(handles.topology_type(1)=='t')
    if(handles.topology_type(8)=='1')
    handles.x1='T1_';
    end
    if(handles.topology_type(8)=='2')
    handles.x1='T2_';
    end
        
end
if(handles.topology_type(1)=='a')
    handles.x1='A_';
end



if (handles.methodSolution(1)=='e') % set elliptic and chebyshev2 as zeross topology
    handles.topology_type = 'zeross';
elseif (handles.methodSolution(1)=='c')
    if (handles.methodSolution(6)=='2')
        handles.topology_type = 'zeross';      
    end
end

if (handles.topology_type(1) == 'z')
    output_img=imread("Z.png");
else
    if(handles.x2=='B')
        %if(((get(handles.stages_pop,'Value') <= ceil(handles.N/2)) && ~mod(N,2)) || ((get(handles.stages_pop,'Value') <= floor(handles.N/2)) && mod(N,2)))
        
        if((get(handles.stages_pop,'Value') == ceil(handles.N/2)) && (floor(handles.N/2))~=ceil(handles.N/2))
            output_img=imread([handles.x1 'B' '.png']);
        else
            if(get(handles.stages_pop,'Value') >= ceil(handles.N/2))
            output_img=imread([handles.x1 'L' '.png']);
            end
            if(get(handles.stages_pop,'Value') <= ceil(handles.N/2))
            output_img=imread([handles.x1 'H' '.png']);
            end
        end
    else
        output_img=imread([handles.x1 handles.x2 '.png']);
    end
end

imshow(output_img);
















% Hints: contents = cellstr(get(hObject,'String')) returns stages_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stages_pop


% --- Executes during object creation, after setting all properties.
function stages_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stages_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function final_results_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes during object creation, after setting all properties.
function text122_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function text22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text122 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in modes_pop.
function modes_pop_Callback(hObject, eventdata, handles)
% hObject    handle to modes_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns modes_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modes_pop


% --- Executes during object creation, after setting all properties.
function modes_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modes_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in methoud_solution_pop.
function methoud_solution_pop_Callback(hObject, eventdata, handles)
% hObject    handle to methoud_solution_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methoud_solution_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methoud_solution_pop





% --- Executes during object creation, after setting all properties.
function table1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modes_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
 %set(handles.table1,'data',[1:15]);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%error

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in methoud_solution_pop.
function table1_Callback(hObject, eventdata, handles)
% hObject    handle to methoud_solution_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns methoud_solution_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from methoud_solution_pop






% --- Executes during object creation, after setting all properties.
function methoud_solution_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to methoud_solution_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in r_c_val_pop.
function r_c_val_pop_Callback(hObject, eventdata, handles)
% hObject    handle to r_c_val_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns r_c_val_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from r_c_val_pop


% --- Executes during object creation, after setting all properties.
function r_c_val_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to r_c_val_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gain_edit_Callback(hObject, eventdata, handles)
% hObject    handle to gain_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gain_edit as text
%        str2double(get(hObject,'String')) returns contents of gain_edit as a double


% --- Executes during object creation, after setting all properties.
function gain_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gain_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
