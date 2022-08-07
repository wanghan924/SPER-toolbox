function varargout = SPER_toolbox(varargin)
%==========================================================================
%
% SPER toolbox developed by Han Wang (zzwang0924@gmail.com),
% Zienkiewicz Centre for Computational Engineering (ZCCE), UK
%
% Please reference to:
% Wang, H., & Xuan, Y. Spatial Pattern Extraction and Recognition (SPER) 
% toolbox for grid-based datasets and its application to machine learning 
% technique. (in preparation)
%
% Wang, H., & Xuan, Y. (2021). Deep learning of extreme rainfall patterns 
% using enhanced spatial random sampling with pattern recognition.
% Hydroinformatics. AGU.
%
% Please contact Han Wang (zzwang0924@gmail.com) with any issue.
%
% Disclaimer:
% This program (hereafter, software) is designed for instructional, educational and research use only.
% Commercial use is prohibited. The software is provided 'as is' without warranty
% of any kind, either express or implied. The software could include technical or other mistakes,
% inaccuracies or typographical errors. The use of the software is done at your own discretion and
% risk and with agreement that you will be solely responsible for any damage and that the authors
% and their affiliate institutions accept no responsibility for errors or omissions in the software
% or documentation. In no event shall the authors or their affiliate institutions be liable to you or
% any third parties for any special, indirect or consequential damages of any kind, or any damages whatsoever.
%
%==========================================================================

%SPER_TOOLBOX MATLAB code file for SPER_toolbox.fig
%      SPER_TOOLBOX, by itself, creates a new SPER_TOOLBOX or raises the existing
%      singleton*.
%
%      H = SPER_TOOLBOX returns the handle to a new SPER_TOOLBOX or the handle to
%      the existing singleton*.
%
%      SPER_TOOLBOX('Property','Value',...) creates a new SPER_TOOLBOX using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to SPER_toolbox_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      SPER_TOOLBOX('CALLBACK') and SPER_TOOLBOX('CALLBACK',hObject,...) call the
%      local function named CALLBACK in SPER_TOOLBOX.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SPER_toolbox

% Last Modified by GUIDE v2.5 30-Jul-2022 11:16:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SPER_toolbox_OpeningFcn, ...
                   'gui_OutputFcn',  @SPER_toolbox_OutputFcn, ...
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


% --- Executes just before SPER_toolbox is made visible.
function SPER_toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output for SPER_toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SPER_toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SPER_toolbox_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Pattern_extraction.
function Pattern_extraction_Callback(hObject, eventdata, handles)
% hObject    handle to Pattern_extraction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1,'reset');
handles.threshold = str2num(get(handles.edit2,'string'));
handles.res = str2num(get(handles.resolution,'string'));
% Check whether the inputs are enough
if isempty(handles.threshold)
   msgbox({'Please input the threshold'});
   return
end

if isempty(handles.data)
   msgbox({'Please load inputs!'});
   return
end

P = handles.data;

switch get(handles.uibuttongroup1,'SelectedObject')%判断选中按钮
    case handles.radiobutton1
        P(P<handles.threshold) = NaN;handles.coreind = 1;
    case handles.radiobutton2
        P(P<=handles.threshold) = NaN;handles.coreind = 1;
    case handles.radiobutton3
        P(P>handles.threshold) = NaN;handles.coreind = 2;
    case handles.radiobutton4
        P(P>=handles.threshold) = NaN;handles.coreind = 2;
end
ind_nan=find(~isnan(P), 1);
if isempty(ind_nan)
   msgbox({'No ROIs are found! Please change threshold!'});
   return
end

axes(handles.axes1);
m = P;
validData = ~isnan(m); % Map of valid (non-nan) locations of matrix m.
[B,~,NN] = bwboundaries(validData);
% imshow(m, [], 'InitialMagnification', 1000); % Display this matrix.
% colormap(gca,'jet');
% cmap = colormap;
% cmap(1,:) = [1 1 1];
% colormap(cmap);
% axis('on', 'image'); % Display tick marks.
f = imagesc(handles.data);
set(f,'alphadata',~isnan(handles.data));

colormap('jet');grid on;
hold on; % Don't let boundaries blow away image because we want to plot them over the image.
i=1;
for k = 1 : length(B)
   thisBoundary = B{k};
   if(k > NN)
     plot(thisBoundary(:,2), thisBoundary(:,1), 'g:','LineWidth',1.5);
   else
     plot(thisBoundary(:,2), thisBoundary(:,1), 'r:','LineWidth',1.5);
     boundaries{i}= thisBoundary;
     i = i+1;
   end
end

[~,I] = sort(cellfun(@length,boundaries),'descend');
boundaries = boundaries(I);kt=1;
for kk = 1:length(boundaries)
    Map = boundaries{kk};
    x = Map(:,1); y = Map(:,2);
    xx=(fix(min(x)):1:fix(max(x)));
    yy=(fix(min(y)):1:fix(max(y)));
    N = 0;
    for i = 1:length(xx)
        for j = 1:length(yy)
            in=inpolygon(xx(i),yy(j),x,y);
            if double(in) == 1 && ~isnan(P(xx(i),yy(j)))
                N = N+1;
                X(N,1) = xx(i);
                Y(N,1) = yy(j);
            end
        end
    end
    if N>5
        XY{kt} = [X,Y];%ik=0;
        for i = 1:length(XY{kt})
            R_data{kt}(i,1) = handles.data(XY{kt}(i,1),XY{kt}(i,2));
        end
        PP = Pattern_extraction([x,y],length(X),XY{kt},R_data{kt},handles.coreind);
        Feature{kt} = PP;
        %PR = R_data{kt}(:);
        R_mean = mean(R_data{kt}(:));
        R_total = sum(R_data{kt}(:))*handles.res*handles.res;
%         for ii = handles.threshold:5:max(R_data{kt})
%             ik = ik+1;
%             R1_mean(1,ik) = mean(PR(PR>ii));
%             R1_total(1,ik) = sum(PR(PR>ii))*sqrt(handles.res);
%         end
        handles.TT(kt,1:11) = [kt PP.YCentriod*handles.res PP.XCentriod*handles.res PP.Area*handles.res*handles.res PP.MajorAxesAngle PP.sp R_mean R_total PP.core(1)*handles.res PP.core(2)*handles.res PP.core(3)];
        handles.Parad{kt,1} = PP.cidx; handles.Parad{kt,2} = PP.indn; handles.Parad{kt,3} = PP.ind1;handles.Parad{kt,4} = PP.cmeans;
        X = [];Y = [];kt=kt+1;%PR = [];
    end
end
handles.ROI_Number = kt-1;
handles.XY = XY;
handles.Feature = Feature;
handles.R_data = R_data;
set(handles.uitable1,'Data',handles.TT);
% Update handles structure
guidata(hObject, handles);clc;

% --- Executes on button press in show.
function show_Callback(hObject, eventdata, handles)
% hObject    handle to show (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes2);
cla(handles.axes2,'reset');
figure
grid on;box on;
XY = handles.XY{handles.ID};
PP = handles.Feature{handles.ID};
RP = handles.R_data{handles.ID};
% Dens = RP / max(RP);
ind = round(RP);
% [~,ind] = sort(Dens);
cmap = jet(length(min(ind):1:max(ind)));
for i = 1:size(XY,1)
    rectangle('Position',[XY(i,2)-0.5,XY(i,1)-0.5,1,1],'FaceColor',cmap(ind(i)-min(ind)+1,1:3),'EdgeColor',cmap(ind(i)-min(ind)+1,1:3));
    hold on
end

xcg=PP.XCentriod;
ycg=PP.YCentriod;
s=PP.s;
xx=[min(PP.boundary(:,1)) max(PP.boundary(:,1))];
yy=s*(xx-xcg)+ycg;
% plot principal axes
plot(yy,xx,'--m','linewidth',1.5);
yy1=[min(PP.boundary(:,2)) max(PP.boundary(:,2))];
xx1 = s*(ycg-yy1)+xcg;
hold on
plot(yy1,xx1,'--m','linewidth',1.5);
hold on
plot(PP.YCentriod,PP.XCentriod,'ok','markerfacecolor','y','markersize',7);
cidx = handles.Parad{handles.ID,1}; 
indn = handles.Parad{handles.ID,2};
ind1 = handles.Parad{handles.ID,3};
cmeans = handles.Parad{handles.ID,4};
for i = 1:sum(indn)
    XY1 = [XY(cidx==ind1(i),1),XY(cidx==ind1(i),2)];
    plot(XY(cidx==ind1(i),2),XY(cidx==ind1(i),1),'kx');
    hold on
    plot(cmeans(ind1(i),2),cmeans(ind1(i),1),'ko','MarkerFaceColor','y');
    hold on
    plot(cmeans(ind1(i),2),cmeans(ind1(i),1),'kx');
    hold on
    % open it if you need calculate the mean and distance between core centre and geometric centre 
%     for j = 1:size(XY1,1)
%         R_data{i}(j,1) = handles.data(XY1(j,1),XY1(j,2));
%     end
%     R_mean{i} = mean(R_data{i}(:));
%     dis{i} = sqrt((PP.YCentriod-cmeans(ind1(i),2)).^2+(PP.XCentriod-cmeans(ind1(i),1)).^2);
    XY1 = [];
end
view(0,90);

axis equal
% camroll(-90)
set(gca,'YDir','reverse','YAxisLocation','right');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName] = uiputfile({'*.jpg','JPEG(*.jpg)';...
                                 '*.bmp','Bitmap(*.bmp)';...
                                 '*.gif','GIF(*.gif)';...
                                 '*.*',  'All Files (*.*)'},...
                                 'Save Picture','Untitled');
PWD1 = pwd;
if FileName==0
   
else
    F = getframe(handles.axes2);
    hf = frame2im(F);
    cd(PathName);
    imwrite(hf, sprintf(FileName));
end
cd(PWD1)
[Filename1,Pathname1]=uiputfile({'*.txt'},'Save Table');
if (Filename1~=0)
  filepath=[Pathname1,'\',Filename1];
else
   return
end
TT=get(handles.uitable1,'Data');

switch get(handles.uibuttongroup1,'SelectedObject')%判断选中按钮
    case handles.radiobutton1
        str = '>=';
    case handles.radiobutton2
        str = '>';
    case handles.radiobutton3
        str = '<=';
    case handles.radiobutton4
        str = '<';
end

fileID=fopen(filepath,'w+');
fprintf(fileID,'%80s \r\n','================================================================================');
fprintf(fileID,'%22s \r\n',strcat('ROIs are extracted ',str,num2str(handles.threshold)));
fprintf(fileID,'%80s \r\n','================================================================================');
fprintf(fileID,'%51s \r\n\r\n','Pattern features are listed below:');
fprintf(fileID,'%6s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s\r\n','No','Location (x)','Location (y)','Size','Orientation(w)','Shape(sp)','Areal average','Total volume','Core(x)','Core(y)','Core(areal average)');
fprintf(fileID,'%6.0f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f %12.2f\r\n',TT');
fprintf(fileID,'%70s \r\n\r\n','');
fclose(fileID);

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
% Select which region to present
choice = get(hObject,'Value');

% Select what copulas to run: based on main functions
if any(choice > handles.ROI_Number)
   msgbox({'Only ',num2str(handles.ROI_Number),' regions have been detected.'...,
       'Please rechoose the region!'});
   return
else
    handles.ID = choice;
end

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Load_Data.
function Load_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.axes1,'reset');
cla(handles.axes2,'reset');
handles.data=[];
[FileName,FilePath] = uigetfile({'*.txt';'*.xlsx';'*.mat'},...
    'Select grid-based inputs (*.txt,*.xlsx,*.mat)','Multiselect','off');
ExPath = fullfile(FilePath, FileName);


if isequal(FileName,0)
   msgbox({'Please load inputs!'});
   return
end

if strcmp(ExPath(end-2:end),'txt') || strcmp(ExPath(end-2:end),'mat')
    % Load data
    handles.data = load(ExPath);
elseif strcmp(ExPath(end-3:end),'xlsx')
    % Load data
    handles.data = xlsread(ExPath);
elseif strcmp(ExPath(end-2:end),'shp')
    handles.data = shaperead(ExPath);
end
if size(handles.data,2)<=2
    msgbox({'The data must have 2 or more than 2 columns'});
    return
end
if ~isempty(handles.data)
    msgbox({'The data have been loaded successfully!'});
else
    msgbox({'Please load inputs!'});
end
% Update handles structure
guidata(hObject, handles);


function resolution_Callback(hObject, eventdata, handles)
% hObject    handle to resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of resolution as text
%        str2double(get(hObject,'String')) returns contents of resolution as a double


% --- Executes during object creation, after setting all properties.
function resolution_CreateFcn(hObject, eventdata, handles)
% hObject    handle to resolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3


% --- Executes during object creation, after setting all properties.
function uibuttongroup1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
