function varargout = InitProcessImageGUI_revised(varargin)
% INITPROCESSIMAGEGUI_REVISED M-file for InitProcessImageGUI_revised.fig
%      INITPROCESSIMAGEGUI_REVISED, by itself, creates a new INITPROCESSIMAGEGUI_REVISED or raises the existing
%      singleton*.
%
%      H = INITPROCESSIMAGEGUI_REVISED returns the handle to a new INITPROCESSIMAGEGUI_REVISED or the handle to
%      the existing singleton*.

%
%      INITPROCESSIMAGEGUI_REVISED('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INITPROCESSIMAGEGUI_REVISED.M with the given input arguments.
%
%      INITPROCESSIMAGEGUI_REVISED('Property','Value',...) creates a new INITPROCESSIMAGEGUI_REVISED or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before InitProcessImageGUI_revised_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to InitProcessImageGUI_revised_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help InitProcessImageGUI_revised

% Last Modified by GUIDE v2.5 20-Sep-2016 09:46:11

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @InitProcessImageGUI_revised_OpeningFcn, ...
    'gui_OutputFcn',  @InitProcessImageGUI_revised_OutputFcn, ...
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


% --- Executes just before InitProcessImageGUI_revised is made visible.
function InitProcessImageGUI_revised_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to InitProcessImageGUI_revised (see VARARGIN)

clear global

% Choose default command line output for InitProcessImageGUI_revised
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes InitProcessImageGUI_revised wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = InitProcessImageGUI_revised_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Process.
function Process_Callback(hObject, eventdata, handles)
% hObject    handle to Process (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --- Executes on selection change in FileList.
function FileList_Callback(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns FileList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileList


% --- Executes during object creation, after setting all properties.
function FileList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Throughput.
function Throughput_Callback(hObject, eventdata, handles)
% hObject    handle to Throughput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in FiberBasisSet.
function FiberBasisSet_Callback(hObject, eventdata, handles)
% hObject    handle to FiberBasisSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat', 'Fiber Basis Set File');
if strcmp(num2str(filename),'0')==0
    set(hObject,'string',[pathname filename]);
else
    set(hObject,'string','-');
end


% --- Executes on button press in CalibrationFile.
function CalibrationFile_Callback(hObject, eventdata, handles)
% hObject    handle to CalibrationFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat', 'Wavelength/Wavenumber Calibration File');
if strcmp(num2str(filename),'0')==0
    set(hObject,'string',[pathname filename]);
else
    set(hObject,'string','-');
end


% --- Executes on button press in initprocessdirectory.
function initprocessdirectory_Callback(hObject, eventdata, handles)
% hObject    handle to initprocessdirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.initprocessdirectory,'string',uigetdir('/Users/chishu/Desktop/data processing', 'Select File Directory'))

% --- Executes on selection change in FileList.
function listbox17_Callback(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FileList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FileList


% --- Executes during object creation, after setting all properties.
function listbox17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FileList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in AddSpec.
function pushbutton97_Callback(hObject, eventdata, handles)
% hObject    handle to AddSpec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in initialprocess.
function initialprocess_Callback(hObject, eventdata, handles)
% hObject    handle to initialprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Options

% Darkspec Filtering
window = -6:6; winsig = 10;

% CCD Image Size
px = 1024; py = 256;%Keren change it to 255 for microscope system

% Aberration Correction
polyorderaberration = 2;
xpixeldrift = 6;
ypixeldrift = 2;

% Throughput Properties
fibernum = 20; % Number of fibers to use in fit - using smaller amount due to only 6 mm (?)
thpeakstripwindow = 3;
thedgedist = 4;

% Neon Properties
npeakstripwindow = 4;
npeaknum = 22;  % this value has to match the length of npeaklambda
npeaklambda = [849.54 859.13 NaN 865.44 NaN 878.06 885.39 886.55 891.95 914.87 920.18 NaN 930.09 NaN 932.65 NaN 942.54 NaN NaN 953.42 954.74 966.54].';

polyorderneon = 3;

% Tylenol Properties
typeakstripwindow = 15;

typeaknum = 21;
typeakwavenum = [NaN 329.2 390.9 465.1 NaN NaN 651.6 710.8 797.2 NaN 857.9 NaN 1168.5 1236.8 NaN 1278.5 1329.9 1371.5 1561.6 NaN 1648.4].';
tyedgedist = 21;

%% Get Files to Process
set(handles.initprocessstatus,'string','Status: Initializing...'); pause(1E-6)
filedir = get(handles.FileDirectory,'string');
s = what(filedir); allfiles = s.mat;
specind = logical(1 - ((1-cellfun('isempty', regexp(allfiles,'throughput'))) + ...
    (1-cellfun('isempty', regexp(allfiles,'neon'))) + ...
    (1-cellfun('isempty', regexp(allfiles,'tylenol'))) + ...
    (1-cellfun('isempty', regexp(allfiles,'darkspec')))));
list = allfiles(specind);

%% Load Dark Spectrum
set(handles.initprocessstatus,'string','Status: Calculating Dark Spectrum...'); pause(1E-6)
darkspec = load([filedir '/darkspec_cali.mat']);
darkspec = squeeze(double(permute(reshape(darkspec.RawData.Spectrum.',px,py, ...
    str2double(darkspec.RawData.NumofKin)),[3 2 1])./...
    str2double(darkspec.RawData.NumofAcu)));
darkspecmad = squeeze(mad(darkspec,1));
darkspecmed = squeeze(median(darkspec,1));
darkspec(darkspec > permute(repmat(darkspecmed+5.*darkspecmad,[1 1 size(darkspec,1)]),[3 1 2])) = nan;
darkspec(darkspec < permute(repmat(darkspecmed-5.*darkspecmad,[1 1 size(darkspec,1)]),[3 1 2])) = nan;
darkspec = squeeze(nanmean(darkspec));

% Remove bad pixels%%%Keren: we do not have bad pixels on new CCD
% 07/14/2020
% darkspec(82:256,185) = mean(darkspec(82:256,[184 186]),2);
% darkspec(115:256,308) = mean(darkspec(115:256,[307 309]),2);
% darkspec(113,501) = mean(darkspec(113,[500 502]),2);

% Filter Dark Spectrum
gausskernel = exp(-(window.^2)./(2.*winsig.^2));
gausskernel = gausskernel./sum(gausskernel(:));
darkspec2 = conv2(darkspec,gausskernel,'same');
mxwindow = max(window);
darkspec2(1:mxwindow,:) = darkspec(1:mxwindow,:); darkspec2((py-mxwindow+1):py,:) = darkspec((py-mxwindow+1):py,:);
darkspec2(:,1:mxwindow) = darkspec(:,1:mxwindow); darkspec2(:,(px-mxwindow+1):px) = darkspec(:,(px-mxwindow+1):px);
% figure; imagesc(darkspec); ct = caxis;
%darkspec = darkspec2; clear darkspec2; %CM comment 12/11/2021
% figure; imagesc(darkspec); caxis(ct);

%% Determine Wavelength Calibration
set(handles.initprocessstatus,'string','Status: Determining Wavelength Calibration...'); pause(1E-6)

neon = load([filedir '/neon.mat']);
accum = str2double(neon.RawData.NumofAcu);
neon = permute(reshape(neon.RawData.Spectrum.',px,py,str2double(neon.RawData.NumofKin)),[2 1 3]);
neon = double(neon);
neon = median(neon,3);
neon = neon-darkspec.*accum;

% % Remove bad pixels %We do not have bad pixels in new CCD %Keren
% neon(82:256,185) = mean(neon(82:256,[184 186]),2);
% neon(115:256,308) = mean(neon(115:256,[307 309]),2);
% neon(113,501) = mean(neon(113,[500 502]),2);

[npeakpixels,npeakheight] = PeakLocationsJ([],sum(neon,1).',npeaknum,npeakstripwindow);

npeaklambdaold = zeros(size(npeaklambda));
while(sum((npeaklambdaold==npeaklambda) + floor((isnan(npeaklambdaold)+isnan(npeaklambda))./2))<npeaknum)
    npeaklambdaold = npeaklambda;
    
    axes(handles.currentstep); cla; plot(sum(neon,1).'); hold on;
    text(npeakpixels,npeakheight,num2str(npeaklambda),'rotation',90);
    axis tight; set(gca,'ytick',[]); box on;
    ylabel('intensity / a.u.'); xlabel('pixel number');
    
    axes(handles.previousstep); cla;
    % for future: make a better way to have this work cross-platform
    figure_handle = openfig('neon.fig'); % figure is opened
    axes_handle = findobj(figure_handle, 'Type', 'Axes'); % find handle to axes in figure
    axes_children_handle = get(axes_handle, 'Children');
    copyobj(axes_children_handle, handles.previousstep); % children of original axes are copied to new axes
    close(figure_handle);
    axes(handles.previousstep);
    axis tight; set(gca,'ytick',[]); box on;
    ylabel('intensity / a.u.'); xlabel('pixel number');
    
    prompt={'Enter the peak wavelength values'};
     
    name='Neon Spectrum';
    numlines=1;
    defaultanswer={num2str(npeaklambda.')};
    
     cellans = inputdlg(prompt,name,numlines,defaultanswer);
     npeaklambda=str2num(cellans{1}).';    
    % Dwight Fairchild
    % This while loop reprompts if the input matrix dimensions do not match
    % the sample's matrix dimensions
    while( ~isequal(size(npeaklambda) ,size(npeaklambdaold)));
        
       prompt={'Enter the peak wavelength values'};
       name='Neon Spectrum';
       numlines=1;
       defaultanswer={num2str(npeaklambda.')};
       
       cellans = inputdlg(prompt,name,numlines,defaultanswer);
       npeaklambda=str2num(cellans{1}).';  
    end
   
   
end
npeaklambda = npeaklambdaold;

ind = logical(1-isnan(npeaklambda));
npeakpixels = npeakpixels(ind); npeaklambda = npeaklambda(ind); npeakheight = npeakheight(ind);

process.wavelength = polyval(polyfit(npeakpixels,npeaklambda,polyorderneon),1:size(neon,2)).';

%% Load Throughput
throughput = load([filedir '/throughput.mat']);
throughput = double(permute(reshape(throughput.RawData.Spectrum.',px,py, ...
    str2double(throughput.RawData.NumofKin)),[3 2 1])./...
    str2double(throughput.RawData.NumofAcu));
if size(throughput,1) > 1
    throughputmad = squeeze(mad(throughput,1));
else
    throughputmad = squeeze(zeros(size(throughput)));
end
throughputmed = squeeze(median(throughput,1));
throughput(throughput > permute(repmat(throughputmed+10.*throughputmad,[1 1 size(throughput,1)]),[3 1 2])) = nan;
throughput(throughput < permute(repmat(throughputmed-10.*throughputmad,[1 1 size(throughput,1)]),[3 1 2])) = nan;
if size(throughput,1) > 1
    throughput = squeeze(nanmean(throughput));
else
    throughput = squeeze(throughput);
end
throughput = throughput-darkspec;

% % Remove bad pixels  %We do not have bad pixels in new CCD %Keren
% throughput(82:256,185) = mean(throughput(82:256,[184 186]),2);
% throughput(115:256,308) = mean(throughput(115:256,[307 309]),2);
% throughput(113,501) = mean(throughput(113,[500 502]),2);

% See NIST ISRM 2241 certificate for polynomial coefficients
A = [1.15596E-17 -9.77171E-14 2.16023E-10 -5.86762E-8 2.28325E-4 9.71937E-2]; % A5 A4 A3 A2 A1 A0
wavenumtemp = (1./785-1./process.wavelength).*10.^7;
ISRM = polyval(A,wavenumtemp);
throughput = throughput./repmat(ISRM.',256,1);%change from 256 to 255  Keren


%% Load Whitelamp 
whitelamp = load([filedir '/whitelamp.mat']);
whitelamp = double(permute(reshape(whitelamp.RawData.Spectrum.',px,py, ...
    str2double(whitelamp.RawData.NumofKin)),[3 2 1])./...
    str2double(whitelamp.RawData.NumofAcu));
if size(whitelamp,1) > 1
    whitelampmad = squeeze(mad(whitelamp,1));
else
    whitelampmad = squeeze(zeros(size(whitelamp)));
end
whitelampmed = squeeze(median(whitelamp,1));
whitelamp(whitelamp > permute(repmat(whitelampmed+10.*whitelampmad,[1 1 size(whitelamp,1)]),[3 1 2])) = nan;
whitelamp(whitelamp < permute(repmat(whitelampmed-10.*whitelampmad,[1 1 size(whitelamp,1)]),[3 1 2])) = nan;
if size(whitelamp,1) > 1
    whitelamp = squeeze(nanmean(whitelamp));
else
    whitelamp = squeeze(whitelamp);
end
whitelamp = whitelamp-darkspec;

%% Determine Wavenumber Calibration
set(handles.initprocessstatus,'string','Status: Determining Wavenumber Calibration...'); pause(1E-6)
tylenol = load([filedir '/tylenol.mat']);
accum = str2double(tylenol.RawData.NumofAcu);
tylenol = permute(reshape(tylenol.RawData.Spectrum.',px,py,str2double(tylenol.RawData.NumofKin)),[2 1 3]);
tylenol = double(tylenol);
tylenol = median(tylenol,3);
tylenol = tylenol-darkspec.*accum;

% % Remove bad pixels  %We do not have bad pixels in new CCD %Keren
% tylenol(82:256,185) = mean(tylenol(82:256,[184 186]),2);
% tylenol(115:256,308) = mean(tylenol(115:256,[307 309]),2);
% tylenol(113,501) = mean(tylenol(113,[500 502]),2);

tylenol = tylenol./repmat(ISRM.',256,1);%Keren from 256
tylenol = double(tylenol)./(repmat(sum(throughput,1),256,1));
stylenol = sum(tylenol,1).';
pixels = (1:px).'; pixels = (pixels-mean(pixels))./std(pixels);
stylenol = stylenol-polyval(polyfit(pixels,stylenol,5),pixels);

[typeakwavelength,typeakheight] = PeakLocationsJ(process.wavelength,stylenol,typeaknum,typeakstripwindow,tyedgedist,0,2);
size(typeakwavenum)
typeakwavenumold= zeros(size(typeakwavenum));
while(sum((typeakwavenumold==typeakwavenum) + floor((isnan(typeakwavenumold)+isnan(typeakwavenum))./2))<typeaknum)
    typeakwavenumold = typeakwavenum;
    axes(handles.currentstep); cla; plot(process.wavelength,stylenol); hold on;
    text(typeakwavelength,typeakheight,num2str(typeakwavenum),'rotation',90);
    axis tight; set(gca,'ytick',[]); box on; xt = xlim;
    ylabel('intensity / a.u.'); xlabel('wavelength / nm');
    
    axes(handles.previousstep); cla;
    % flag to fix this also
    figure_handle = openfig('tylenol.fig'); % figure is opened
    axes_handle = findobj(figure_handle, 'Type', 'Axes'); % find handle to axes in figure
    axes_children_handle = get(axes_handle, 'Children');
    copyobj(axes_children_handle, handles.previousstep); % children of original axes are copied to new axes
    close(figure_handle);
    axes(handles.previousstep);
    axis tight; set(gca,'ytick',[]); box on; xlim(xt);
    ylabel('intensity / a.u.'); xlabel('wavelength / nm');
    
    prompt = {'Enter the peak wavenumber values'};
    name = 'Tylenol Spectrum';
    numlines =1 ;
    defaultanswer = {num2str(typeakwavenum.')};
    cellans = inputdlg(prompt,name,numlines,defaultanswer);
    typeakwavenum = str2num(cellans{1}).';
    
    % Dwight Fairchild
    % This while loop reprompts if the input matrix dimensions do not match
    % the sample's matrix dimensions
     while( ~isequal(size(typeakwavenum) ,size(typeakwavenumold)));
        
        prompt={'Enter the peak wavelength values'};
        name='Neon Spectrum';
        numlines=1;
        defaultanswer={num2str(typeakwavenum.')};
       
      cellans = inputdlg(prompt,name,numlines,defaultanswer);
       typeakwavenum=str2num(cellans{1}).';  
    end
end

ind = logical(1-isnan(typeakwavenum));
typeakwavelength = typeakwavelength(ind); typeakwavenum = typeakwavenum(ind);

process.laserwavelength = mean(((typeakwavenum*10^-7) + 1./typeakwavelength).^-1);
process.wavenum = (1./process.laserwavelength - 1./process.wavelength)*10^7;

%% Determine Aberration Correction
set(handles.initprocessstatus,'string','Status: Determining Aberration Correction...'); pause(1E-6)
% Correction along CCD height
[idealycps,peakheight] = PeakLocationsJ([],sum(throughput,2),fibernum,thpeakstripwindow,thedgedist);
ridealycps = round(idealycps);

for ijk = 1:fibernum
    [a,measuredycps(ijk,:)] = max(throughput( (-ypixeldrift:ypixeldrift) + ridealycps(ijk),:));
end
measuredycps = measuredycps+repmat(ridealycps,1,px)-(ypixeldrift+1);

pixels = (1:px).'; pixels = (pixels-mean(pixels))./std(pixels);
% Initialize polynomial matrix
polynom = zeros(length(pixels),polyorderaberration+1);
% Values range from 1 to 2^order
for j=0:1:polyorderaberration
    polynom(:,j+1)=pixels.^j;
end

% figure; imagesc(throughput); colormap('gray'); hold on; plot(measuredycps.','r','linewidth',1)
[a,measuredycps] = OLSJ(measuredycps.',polynom); measuredycps = measuredycps.'; % replace with polyval and polyfit
% figure; imagesc(throughput); colormap('gray'); hold on; plot(measuredycps.','r','linewidth',1)

for ijk = 1:px
    Yi(:,ijk) = polyval(polyfit(idealycps,measuredycps(:,ijk),polyorderaberration),1:py);
end

idealxcps = PeakLocationsJ([],stylenol,typeaknum,typeakstripwindow,tyedgedist,0,2);
ridealxcps = round(idealxcps);

for ijk = 1:typeaknum
    [a,measuredxcps(ijk,:)] = max(tylenol(:, (-xpixeldrift:xpixeldrift) + ridealxcps(ijk)),[],2);
end
measuredxcps = measuredxcps+repmat(ridealxcps,1,py)-(xpixeldrift+1);

pixels = (1:py).'; pixels = (pixels-mean(pixels))./std(pixels);
% Initialize polynomial matrix
polynom = zeros(length(pixels),polyorderaberration+1);
% Values range from 1 to 2^order
for j=0:1:polyorderaberration
    polynom(:,j+1)=pixels.^j;
end
% figure; imagesc(tylenol); colormap('gray'); hold on; plot(measuredxcps.',1:255,'r','linewidth',1)
[a,measuredxcps] = OLSJ(measuredxcps.',polynom); measuredxcps = measuredxcps.'; % replace with polyval and polyfit
% figure; imagesc(tylenol); colormap('gray'); hold on; plot(measuredxcps.',1:256,'r','linewidth',1)

for ijk = 1:py
    Xi(ijk,:) = polyval(polyfit(idealxcps,measuredxcps(:,ijk),polyorderaberration),1:px);
end

%% Determine crop pixels
ind = sum(Xi>=1); mx = max(ind);
cpx1 = find(ind==mx,1,'first');
ind = sum(Xi<=1024); mx = max(ind);
cpx2 = find(ind==mx,1,'last');
ind = sum(Yi>=1,2); mx = max(ind);
cpy1 = find(ind==mx,1,'first');
ind = sum(Yi<=256,2); mx = max(ind);
cpy2 = find(ind==mx,1,'last');
cpx = cpx2-cpx1+1;
cpy = cpy2-cpy1+1;
process.wavenum = process.wavenum(cpx1:cpx2);
process.wavelength = process.wavelength(cpx1:cpx2);

%% Determine Fiber Basis   %No longer doing this CM 12/11/2021
% set(handles.initprocessstatus,'string','Status: Determining Fiber Basis...'); pause(1E-6)
% throughput = throughput.*repmat(ISRM.',256,1); % Return to uncorrected throughput
% ISRM = ISRM(cpx1:cpx2);
% throughput = interp2(throughput,Xi,Yi,'spline'); % Apply aberration correction
% throughput = throughput(cpy1:cpy2,cpx1:cpx2)./repmat(ISRM.',cpy,1); % Crop spectrum and correct for throughput
% sthroughput = sum(throughput,2);
% [peaklocs,peakheight] = PeakLocationsJ([],sthroughput,fibernum,thpeakstripwindow,thedgedist-py+cpy2);
% [~,~,~,~,fiberbasistot] = GFitFiberBasisJ((1:cpy).',sthroughput,fibernum,peaklocs,peakheight);
% 
% % TolX = 1e-16; % 1e-6 by default
% % options = optimset('TolX',TolX);
% for ijk = 1:size(throughput,2)
%     %     c = lsqnonneg(fiberbasistot,throughput(:,ijk),options);
%     %     fiberbasis{ijk} = fiberbasistot.*repmat(c.',cpy,1);
%     [~,~,fiberbasis(ijk)] = OLSJ(throughput(:,ijk),fiberbasistot);
% end
% for kkk=1:40
%     fiberbasis2(:,kkk)=fiberbasis{1}(:,kkk)./max(fiberbasis{1}(:,kkk)).*max(fiberbasis{1}(:,1));
% end
% 
%% Clear unused variables
%clearvars -except hObject eventdata handles process filedir fiberbasis2 list darkspec Xi Yi fiberbasis px py cpx1 cpx2 cpy1 cpy2 cpx cpy
clearvars -except hObject eventdata handles process filedir list darkspec Xi Yi px py cpx1 cpx2 cpy1 cpy2 cpx cpy throughput ISRM whitelamp

%% Process Data
window = -6:6; winsig = 10;
set(handles.initprocessstatus,'string','Status: Calculating Dark Spectrum...'); pause(1E-6)
darkspec = load([filedir '/darkspec.mat']);
darkspec = squeeze(double(permute(reshape(darkspec.RawData.Spectrum.',px,py, ...
    str2double(darkspec.RawData.NumofKin)),[3 2 1])./...
    str2double(darkspec.RawData.NumofAcu)));
darkspecmad = squeeze(mad(darkspec,1));
darkspecmed = squeeze(median(darkspec,1));
darkspec(darkspec > permute(repmat(darkspecmed+5.*darkspecmad,[1 1 size(darkspec,1)]),[3 1 2])) = nan;
darkspec(darkspec < permute(repmat(darkspecmed-5.*darkspecmad,[1 1 size(darkspec,1)]),[3 1 2])) = nan;
darkspec = squeeze(nanmean(darkspec));

% % Remove bad pixels %We do not have bad pixels in new CCD %Keren
% darkspec(82:256,185) = mean(darkspec(82:256,[184 186]),2);
% darkspec(115:256,308) = mean(darkspec(115:256,[307 309]),2);
% darkspec(113,501) = mean(darkspec(113,[500 502]),2);

% Filter Dark Spectrum
gausskernel = exp(-(window.^2)./(2.*winsig.^2));
gausskernel = gausskernel./sum(gausskernel(:));
darkspec2 = conv2(darkspec,gausskernel,'same');
mxwindow = max(window);
darkspec2(1:mxwindow,:) = darkspec(1:mxwindow,:); darkspec2((py-mxwindow+1):py,:) = darkspec((py-mxwindow+1):py,:);
darkspec2(:,1:mxwindow) = darkspec(:,1:mxwindow); darkspec2(:,(px-mxwindow+1):px) = darkspec(:,(px-mxwindow+1):px);
% figure; imagesc(darkspec); ct = caxis;
%darkspec = darkspec2; clear darkspec2; %CM 12/11/2021

%-----------Ab. correct and crop calibration samples white lamp and green glass  %CM
%12/11/2021
whitelamp = interp2(whitelamp,Xi,Yi,'spline'); 
whitelamp = whitelamp(cpy1:cpy2,cpx1:cpx2); %cpx1 => 10 CM 12/12/2021

throughput = interp2(throughput,Xi,Yi,'spline'); 
throughput = throughput(cpy1:cpy2,cpx1:cpx2); %cpx1 => 10 CM 12/12/2021
ISRM2 = ISRM(cpx1:cpx2); 



set(handles.initprocessstatus,'string','Status: Processing Data...'); pause(1E-6)
if isempty(list) == 0
    if isdir([filedir '/initprocess/']) == 0
        mkdir(filedir,'initprocess');
    end
    process.location = []; process.framenum = []; process.exptime = [];
    process.spec = []; process.shotnoise = []; process.accumnum = [];
    process.saturated = [];
    process = orderfields(process);
    
    totiter = length(list);
    curriter = 1;
    tic
    for ijk = 1:length(list) % For each spectral data file
        savefilename = [filedir '/initprocess/' list{ijk}];
        if exist(savefilename,'file')~=2
            % read data file
            process.location = [filedir '/' list{ijk}];
            load(process.location);
            
            process.framenum = str2double(RawData.NumofKin);
            process.exptime = str2double(RawData.ExposureTime);
            process.accumnum = str2double(RawData.NumofAcu);
            process.spec = [];
            process.shotnoise = [];
            
            if sum(RawData.Spectrum(:) == 65535) > 60% if more than one saturated pixel
                process.saturated = 1;
                disp([process.location ' image saturated!']);
            else
                process.saturated = 0;
                Z = permute(reshape(RawData.Spectrum.',px,py,str2double(RawData.NumofKin)),[2 1 3]);%keren just a mark, I forgot why
                clear RawData
                Z = double(Z);
                Z = Z-repmat(darkspec,[1 1 size(Z,3)]).*process.accumnum; % Subtract dark spectrum
                Z = Z.*7; % Convert to photoelectrons
                
                
%-----------manual set the rows for each leg CM 12/11/2021 not doing
%fiberbasis to get fiber spectrum
                leg{1} = [67:95];  %0mm offset (4 fibers)              
                leg{2} = [1:66];   %3mm offset (12 fibers)
                leg{3} = [96:252]; %6mm offset (26 fibers)
%-----------manual set the rows for each leg CM 12/11/2021

                for klm = 1:size(Z,3)
%                     % Remove bad pixels  %Keren
%                     Z(82:256,185,klm) = mean(Z(82:256,[184 186],klm),2);
%                     Z(115:256,308,klm) = mean(Z(115:256,[307 309],klm),2);
%                     Z(113,501,klm) = mean(Z(113,[500 502],klm),2);
                    
                    % Aberration Correction
                    ZCt = interp2(Z(:,:,klm),Xi,Yi,'spline');
                    
                    % Crop image
                    ZC = ZCt(cpy1:cpy2,cpx1:cpx2); 
                    process.image(:,:,klm)=ZC;
                    
                    
                    for i = 1:size(leg,2)
                       
                        throughput_leg(:,i) = sum(throughput(leg{i},:)); %get throughput leg CM 12/11/2021
                        throughput_leg_s(:,i) = smooth(throughput_leg(:,i),100); %smooth to calculate spectral response CM 12/11/2021
                        throughput_leg_s(:,i) = throughput_leg_s(:,i)./ISRM2; %calculate spectral response from smooth green glass and NIST standard CM 12/11/2021
                        
                        ZC_leg(:,i) = sum(ZC(leg{i},:)); %get data leg CM 12/11/2021
                        ZC_leg_thru(:,i) = ZC_leg(:,i)./throughput_leg_s(:,i);  %correct for spectral throughput CM 12/11/2021
                        
                        whitelamp_leg(:,i) = sum(whitelamp(leg{i},:));  %get whitelamp leg 12/11/2021
                        s = smooth(whitelamp_leg(:,i),50); %CM
                        fp = whitelamp_leg(:,i)./s; %low pass filter to get fixed pattern CM 12/11/2021
%                         
                        
                        ZC_2(:,i) = ZC_leg_thru(:,i)./(fp); %
                        %ZC_2(:,i) = ZC_leg_thru(:,i); 
                    end
                    ZC_3(:,:,klm) = ZC_2; 
            
                    
                    %-----No longer extracting fiber spectra CM 12/11/2021
%                     % Extract spectrum from each fiber
%                     for jkl = 1:cpx
%                         [process.spec(jkl,:,klm)] = ...
%                             OLSJ(ZC(:,jkl),[fiberbasis{jkl}]);
%                         
%                         % Covariance matrix of noise (assuming measurement noise is uncorrelated)
%                         % See Determination of uncertainty in parameters extracted
%                         % from single spectroscopic measurements, equation (4)
%                         Cw = diag(ZC(:,jkl));
%                         variance = diag(inv(fiberbasis{jkl}'/Cw*fiberbasis{jkl}));
%                         variance(variance<0) = 0; variance(isnan(variance)) = 0;
%                         
%                         process.shotnoise(jkl,:,klm) = sqrt(variance);
%                         
%                         % Add readout noise?
%                         
%                     end
                    
                    %-------scaling factor from using fiberbasis => no
                    %longer needed I think CM 12/11/2021
%                     normfac = repmat(mean(process.shotnoise(:,:,klm)),size(process.shotnoise,1),1);
%                     process.spec(:,:,klm) = process.spec(:,:,klm)./normfac;%filt shot noise
%                     meansig = repmat(mean(process.spec(:,:,klm)),size(process.spec,1),1);%mean spectrum
%                     process.spec(:,:,klm) = process.spec(:,:,klm).*meansig;
%                     process.shotnoise(:,:,klm) = process.shotnoise(:,:,klm)./normfac.*meansig;
                end
                process.spec = ZC_3;
            end
            
            % Save Data
            save(savefilename,'process')
            
            % Display remaining time
            [hours,minutes,seconds] = RemainingTimeJ(curriter,totiter);
            minutes = minutes+hours.*60;
            set(handles.initprocessstatus,'string',['Status: Processing Data... ' ...
                num2str(minutes) ' min ' num2str(seconds) ' sec remaining']); pause(1E-6)
            curriter = curriter+1;
        else
            totiter = totiter-1;
            tic
        end
    end
    
    set(handles.initprocessstatus,'string','Status: Ready');
else
    errordlg('No Data Files', 'Error');
end


% --- Executes on button press in FiberBasisSet.
function pushbutton84_Callback(hObject, eventdata, handles)
% hObject    handle to FiberBasisSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in CalibrationFile.
function wavenumbercal_Callback(hObject, eventdata, handles)
% hObject    handle to CalibrationFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in svgsmoothonoff.
function svgsmoothonoff_Callback(hObject, eventdata, handles)
% hObject    handle to svgsmoothonoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value') == 1
    set(handles.smoothpolyorder,'enable','on'); set(handles.smoothwindowsize,'enable','on');
else
    set(handles.smoothpolyorder,'enable','off'); set(handles.smoothwindowsize,'enable','off');
end

% Hint: get(hObject,'Value') returns toggle state of svgsmoothonoff


function bkgdpolyorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to polyorderaberration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of polyorderaberration as text
%        str2double(get(hObject,'String')) returns contents of polyorderaberration as a double


function polyorderaberration_Callback(hObject, eventdata, handles)
% hObject    handle to polyorderaberration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of polyorderaberration as text
%        str2double(get(hObject,'String')) returns contents of polyorderaberration as a double


% --- Executes during object creation, after setting all properties.
function polyorderaberration_CreateFcn(hObject, eventdata, handles)
% hObject    handle to polyorderaberration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c1_Callback(hObject, eventdata, handles)
% hObject    handle to c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c1 as text
%        str2double(get(hObject,'String')) returns contents of c1 as a double


% --- Executes during object creation, after setting all properties.
function c1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c2_Callback(hObject, eventdata, handles)
% hObject    handle to c2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c2 as text
%        str2double(get(hObject,'String')) returns contents of c2 as a double


% --- Executes during object creation, after setting all properties.
function c2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c1t_Callback(hObject, eventdata, handles)
% hObject    handle to c1t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c1t as text
%        str2double(get(hObject,'String')) returns contents of c1t as a double


% --- Executes during object creation, after setting all properties.
function c1t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c1t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c2t_Callback(hObject, eventdata, handles)
% hObject    handle to c2t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c2t as text
%        str2double(get(hObject,'String')) returns contents of c2t as a double


% --- Executes during object creation, after setting all properties.
function c2t_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c2t (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function crrmad_Callback(hObject, eventdata, handles)
% hObject    handle to crrmad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of crrmad as text
%        str2double(get(hObject,'String')) returns contents of crrmad as a double


% --- Executes during object creation, after setting all properties.
function crrmad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to crrmad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tempprocess.
function finalprocess_Callback(hObject, eventdata, handles)
% hObject    handle to tempprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

process.filedirectory = get(handles.initprocessdirectory,'string');
s = what(process.filedirectory ); process.list = s.mat;
[process.list] = sort_nat(process.list).';

% Wavenumber regions
c1 = str2double(get(handles.c1,'string')); c2 = str2double(get(handles.c2,'string'));
c1t = str2double(get(handles.c1t,'string')); c2t = str2double(get(handles.c2t,'string'));
c1n = str2double(get(handles.c1n,'string')); c2n = str2double(get(handles.c2n,'string'));
process.processoptions.c1 = c1; process.processoptions.c2 = c2;
process.processoptions.c1t = c1t; process.processoptions.c2t = c2t;
process.wavenum = (process.processoptions.c1:2:process.processoptions.c2).';

normflag = get(handles.normonoff,'value');
if normflag == 1
    process.processoptions.c1n = c1n; process.processoptions.c2n = c2n;
    c1n = find(process.wavenum==c1n); c2n = find(process.wavenum==c2n);
end
c1t = find(process.wavenum==c1t); c2t = find(process.wavenum==c2t);
c1 = find(process.wavenum==c1); c2 = find(process.wavenum==c2);

if strcmp(get(handles.addlb,'string'),'-')==1
    addlb = [];
else
    addlb = load(get(handles.addlb,'string'));
    addlb = addlb.basis;
end

process.processoptions.polyorder = str2double(get(handles.bkgdpolyorder,'string'));
process.processoptions.iter = str2double(get(handles.iter,'string'));
%process.processoptions.iter = 100; %CM 4/19
process.processoptions.crrmad = str2double(get(handles.crrmad,'string'));
process.processoptions.convergepercent = str2double(get(handles.convergepercent,'string'));
process.processoptions.peakremovalflag = get(handles.peakremovalflag,'value');
process.processoptions.noisefac = str2double(get(handles.noisefac,'string'));
process.processoptions.noisecalc = get(handles.noisecalc,'value');

tic
for jkl = 1:size(process.list,2)
    initprocess = load([process.filedirectory '/' process.list{jkl}]);
    process.laserwavelength{jkl} = initprocess.process.laserwavelength;
    process.wavelength{jkl} = (1./process.laserwavelength{jkl} - process.wavenum.*1E-7).^-1;
    %process.image{jkl}=initprocess.process.image;Keren image
    
    if isfield(initprocess.process,'saturated') == 1
        process.saturated(jkl) = initprocess.process.saturated;
    else
        process.saturated(jkl) = 0;
    end
    
    if process.saturated(jkl) == 0
        spec = interp1(initprocess.process.wavelength,initprocess.process.spec,process.wavelength{jkl},'spline');
        %-CM comment 12/13/2021
        %shotnoise = interp1(initprocess.process.wavelength,initprocess.process.shotnoise,process.wavelength{jkl},'spline');
        
% %-------------no longer using shot noise CM 12/13/2021
%         [process.spec(jkl),process.shotnoise(jkl),process.empiricalnoise(jkl)] = ...
%             CRRPolyFit(process.wavenum,spec,shotnoise ,...
%             process.processoptions.polyorder,c1,c2, ...
%             process.processoptions.crrmad,addlb);
%-------------------------- CM 12/13/2021
        [process.spec(jkl)] = ...
            CRRPolyFit(process.wavenum,spec,...
            process.processoptions.polyorder,c1,c2, ...
            process.processoptions.crrmad,addlb);
        
        
        
        if get(handles.svgsmoothonoff,'value') == 1 % Smoothing
            process.processoptions.SVGsmooth = 'on';
            process.processoptions.SVGsmoothpolyorder = ...
                str2double(get(handles.smoothpolyorder,'string')); % Polynomial Order;
            process.processoptions.SVGsmoothwindowsize = ...
                str2double(get(handles.smoothwindowsize,'string')); % Window Size;
            for nop = 1:size(process.spec{jkl},2)
                process.spec{jkl}(:,nop) = SvgSmoothJ( ...
                    process.processoptions.SVGsmoothpolyorder, ...
                    process.processoptions.SVGsmoothwindowsize, ...
                    process.spec{jkl}(:,nop));
                    
            end
        else
            process.processoptions.SVGsmooth = 'off';
            process.processoptions.SVGsmoothpolyorder = [];
            process.processoptions.SVGsmoothwindowsize = [];
        end
    else
        process.spec{jkl} = nan([length(process.wavelength{jkl}),40]);
        process.shotnoise{jkl} = process.spec{jkl}; %unsure CM 12/13
        process.empiricalnoise{jkl} = process.spec{jkl}; %unsure CM 12/13
    end
    
    %-not needed already have legs for offsets CM 12/13/2021
%     [process.specradii{jkl},process.shotnoiseradii{jkl},process.empiricalnoiseradii{jkl}] = ...
%         ConvertFiberDataRadii(process.spec{jkl},process.shotnoise{jkl},process.empiricalnoise{jkl});
%     
    if process.processoptions.noisecalc == 1
        [process.anitaspec{jkl}] = IanitaJ(process.wavenum, ...
            process.spec{jkl},process.processoptions.polyorder, ...
            c1,c2,process.processoptions.iter,addlb, ...
            process.processoptions.convergepercent, ...
            process.processoptions.peakremovalflag, ...
            process.processoptions.noisefac);
        
        %-not needed CM 12/13/2021
%         [process.anitaspecradii{jkl}] = IanitaJ(process.wavenum, ...
%             process.specradii{jkl},process.processoptions.polyorder, ...
%             c1,c2,process.processoptions.iter,addlb, ...
%             process.processoptions.convergepercent, ...
%             process.processoptions.peakremovalflag, ...
%             process.processoptions.noisefac);
    else
        
        %no shot noise CM 12/13/2021
%         [process.anitaspec{jkl}] = IanitaJ(process.wavenum, ...
%             process.spec{jkl},process.processoptions.polyorder, ...
%             c1,c2,process.processoptions.iter,addlb, ...
%             process.processoptions.convergepercent, ...
%             process.processoptions.peakremovalflag, ...
%             process.processoptions.noisefac, ...
%             process.shotnoise{jkl});
        
            [process.anitaspec{jkl}] = IanitaJ(process.wavenum, ...
            process.spec{jkl},process.processoptions.polyorder, ...
            c1,c2,process.processoptions.iter,addlb, ...
            process.processoptions.convergepercent, ...
            process.processoptions.peakremovalflag, ...
            process.processoptions.noisefac); %CM 12/13/2021
        
        %not needed CM 12/13/2021
%         [process.anitaspecradii{jkl}] = IanitaJ(process.wavenum, ...
%             process.specradii{jkl},process.processoptions.polyorder, ...
%             c1,c2,process.processoptions.iter,addlb, ...
%             process.processoptions.convergepercent, ...
%             process.processoptions.peakremovalflag, ...
%             process.processoptions.noisefac, ...
%             process.shotnoiseradii{jkl});
    end
    
    % Normalization
    if normflag == 1
        process.processoptions.normalization = 'not yet operational';
    else
        process.processoptions.normalization = 'off';
    end
    [hours, minutes, seconds] = RemainingTimeJ(jkl,size(process.list,2));
    minutes = minutes+hours.*60;
    set(handles.finalprocessstatus,'string',['Status: Processing Data... ' ...
        num2str(minutes) ' min ' num2str(seconds) ' sec remaining']); pause(1E-6)
end

% tempprocess = process;
% fields{1} = 'spec'; fields{2} = 'shotnoise'; fields{3} = 'wavenum';
% fields{4} = 'wavelength'; fields{5} = 'list'; fields{6} = 'anitaspec';
% fields{7} = 'specradii'; fields{8} = 'anitaspecradii'; fields{9} = 'empiricalnoise';
% process = rmfield(process,fields);

clearvars -except handles process c1 c2 c1t c2t addlb

% Calculate mean spectrum and shotnoise
process.meanspec = cell2mat(cellfun(@mean,process.spec, ...
    num2cell(ones(size(process.spec)).*2),'uniformoutput',0));

%---------CM 12/13/2021
% %specnum = cell2mat(cellfun(@size,process.shotnoise, ...
%     num2cell(ones(size(process.shotnoise)).*2),'uniformoutput',0));
% 
%     process.meanshotnoise = cellfun(@power,process.shotnoise, ...
%     num2cell(ones(size(process.shotnoise)).*2),'uniformoutput',0);
% process.meanshotnoise = cell2mat(cellfun(@sum,process.meanshotnoise, ...
%     num2cell(ones(size(process.meanshotnoise)).*2),'uniformoutput',0));
% process.meanshotnoise = repmat(1./specnum,size(process.meanshotnoise,1),1).*sqrt(process.meanshotnoise);
% 
% process.meanempiricalnoise = cellfun(@power,process.empiricalnoise, ...
%     num2cell(ones(size(process.empiricalnoise)).*2),'uniformoutput',0);
% process.meanempiricalnoise = cell2mat(cellfun(@sum,process.meanempiricalnoise, ...
%     num2cell(ones(size(process.meanempiricalnoise)).*2),'uniformoutput',0));
% process.meanempiricalnoise = repmat(1./specnum,size(process.meanempiricalnoise,1),1).*sqrt(process.meanempiricalnoise);
%--------CM 12/13/2021

if process.processoptions.noisecalc == 1
    [process.meananitaspec] = IanitaJ(process.wavenum, ...
        process.meanspec,process.processoptions.polyorder, ...
        c1,c2,process.processoptions.iter,addlb, ...
        process.processoptions.convergepercent, ...
        process.processoptions.peakremovalflag, ...
        process.processoptions.noisefac);
else
    %---CM 12/13/2021
%     [process.meananitaspec] = IanitaJ(process.wavenum, ...
%         process.meanspec,process.processoptions.polyorder, ...
%         c1,c2,process.processoptions.iter,addlb, ...
%         process.processoptions.convergepercent, ...
%         process.processoptions.peakremovalflag, ...
%         process.processoptions.noisefac, ...
%         process.meanshotnoise);

        [process.meananitaspec] = IanitaJ(process.wavenum, ...
        process.meanspec,process.processoptions.polyorder, ...
        c1,c2,process.processoptions.iter,addlb, ...
        process.processoptions.convergepercent, ...
        process.processoptions.peakremovalflag, ...
        process.processoptions.noisefac); %CM 12/13/2021
end

% Cropping
for jkl = 1:length(process.list)
    process.spec{jkl} = process.spec{jkl}(c1t:c2t,:);
    %process.specradii{jkl} = process.specradii{jkl}(c1t:c2t,:);   %CM
    %12/13/2021
    
    process.anitaspec{jkl} = process.anitaspec{jkl}(c1t:c2t,:);
    %process.anitaspecradii{jkl} = process.anitaspecradii{jkl}(c1t:c2t,:);
    %CM 12/13/2021
    
    %process.shotnoise{jkl} = process.shotnoise{jkl}(c1t:c2t,:); %CM
    %12/13/2021
    
    %process.empiricalnoise{jkl} = process.empiricalnoise{jkl}(c1t:c2t,:);
    %CM 12/13/2021
    %process.shotnoiseradii{jkl} = process.shotnoiseradii{jkl}(c1t:c2t,:);
    %CM 12/13/2021
    
    %process.empiricalnoiseradii{jkl} =
    %process.empiricalnoiseradii{jkl}(c1t:c2t,:); %CM 12/13/2021
    process.wavelength{jkl} = process.wavelength{jkl}(c1t:c2t);
end
process.meanspec = process.meanspec(c1t:c2t,:);
%process.meanshotnoise = process.meanshotnoise(c1t:c2t,:); %CM 12/13/2021
%process.meanempiricalnoise = process.meanempiricalnoise(c1t:c2t,:); %CM
%12/13/2021
process.meananitaspec = process.meananitaspec(c1t:c2t,:);
process.wavenum = process.wavenum(c1t:c2t);
process.list = process.list;

axes(handles.previousstep); cla;
plotR(process.wavenum,process.meanspec./ ...
    repmat(sum(process.meanspec),size(process.meanspec,1),1));
axes(handles.currentstep); cla;
plotR(process.wavenum,process.meananitaspec./ ...
    repmat(mad(process.meananitaspec),size(process.meananitaspec,1),1));

prompt={'Comment:'};
name='Comment';
defaultanswer={'none'};
comment=inputdlg(prompt,name,1,defaultanswer);
process.comment = comment;
process = orderfields(process);
process.processoptions = orderfields(process.processoptions);
uisave('process')

set(handles.finalprocessstatus,'string','Status: Ready');


function [specout,shotnoiseout,empiricalnoiseout] = ...
    ConvertFiberDataRadii(specin,shotnoisein,empiricalnoisein)

% 5 radii
% % % ind{1} = 19;
% % % ind{2} = [1 2 4 8 16 17];
% % % ind{3} = [3 5 6 7 9 10 11 12 13 14 15 18];
% % % ind{4} = [20 24 27 29 31 32 33 37 38 40];
% % % ind{5} = [21 22 23 25 26 28 30 34 35 36 39];

% 3 radii
ind{1} = [1 2 4 8 16 17 19];
ind{2} = [3 5 6 7 9 10 11 12 13 14 15 18];
ind{3} = 20:40;

% % % 3 radii
% % % ind{1} = 19;
% % % ind{2} = [1:19 20 24 27 29 31 32 33 37 38 40];
% % % ind{3} = [21 22 23 25 26 28 30 34 35 36 39];

% % % 3 radii
% % % ind{1} = 19;
% % % ind{2} = 1:18;
% % % ind{3} = 20:40;

% % % 2 radii
% ind{1} = [1 2 4 8 16 17 19];
% ind{2} = [3 5 6 7 9 10 11 12 13 14 15 18 20:40];

specout = zeros(size(specin,1),size(ind,2));
shotnoiseout = zeros(size(specin,1),size(ind,2));
empiricalnoiseout = zeros(size(specin,1),size(ind,2));

for ijk = 1:size(ind,2)
    normfac = size(specin(:,ind{ijk}),2);
    specout(:,ijk) = sum(specin(:,ind{ijk}),2)./normfac;
    shotnoiseout(:,ijk) = sqrt(sum(shotnoisein(:,ind{ijk}).^2,2))./normfac;
    empiricalnoiseout(:,ijk) = sqrt(sum(empiricalnoisein(:,ind{ijk}).^2,2))./normfac;
end

% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in normonoff.
function normonoff_Callback(hObject, eventdata, handles)
% hObject    handle to normonoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value') == 1
    set(handles.c1n,'enable','on'); set(handles.c2n,'enable','on');
else
    set(handles.c1n,'enable','off'); set(handles.c2n,'enable','off');
end
% Hint: get(hObject,'Value') returns toggle state of normonoff



function c1n_Callback(hObject, eventdata, handles)
% hObject    handle to c1n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c1n as text
%        str2double(get(hObject,'String')) returns contents of c1n as a double


% --- Executes during object creation, after setting all properties.
function c1n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c1n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function c2n_Callback(hObject, eventdata, handles)
% hObject    handle to c2n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of c2n as text
%        str2double(get(hObject,'String')) returns contents of c2n as a double


% --- Executes during object creation, after setting all properties.
function c2n_CreateFcn(hObject, eventdata, handles)
% hObject    handle to c2n (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in removeoutlierspeconoff.
function removeoutlierspeconoff_Callback(hObject, eventdata, handles)
% hObject    handle to removeoutlierspeconoff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of removeoutlierspeconoff



function smoothwindowsize_Callback(hObject, eventdata, handles)
% hObject    handle to smoothwindowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smoothwindowsize as text
%        str2double(get(hObject,'String')) returns contents of smoothwindowsize as a double


% --- Executes during object creation, after setting all properties.
function smoothwindowsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothwindowsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function compilespecregexp_Callback(hObject, eventdata, handles)
% hObject    handle to compilespecregexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of compilespecregexp as text
%        str2double(get(hObject,'String')) returns contents of compilespecregexp as a double


% --- Executes during object creation, after setting all properties.
function compilespecregexp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to compilespecregexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function convergepercent_Callback(hObject, eventdata, handles)
% hObject    handle to convergepercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of convergepercent as text
%        str2double(get(hObject,'String')) returns contents of convergepercent as a double


% --- Executes during object creation, after setting all properties.
function convergepercent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to convergepercent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in peakremovalflag.
function peakremovalflag_Callback(hObject, eventdata, handles)
% hObject    handle to peakremovalflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of peakremovalflag



function iter_Callback(hObject, eventdata, handles)
% hObject    handle to iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iter as text
%        str2double(get(hObject,'String')) returns contents of iter as a double


% --- Executes during object creation, after setting all properties.
function iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function noisefac_Callback(hObject, eventdata, handles)
% hObject    handle to noisefac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of noisefac as text
%        str2double(get(hObject,'String')) returns contents of noisefac as a double


% --- Executes during object creation, after setting all properties.
function noisefac_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noisefac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in noisecalc.
function noisecalc_Callback(hObject, eventdata, handles)
% hObject    handle to noisecalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns noisecalc contents as cell array
%        contents{get(hObject,'Value')} returns selected item from noisecalc


% --- Executes during object creation, after setting all properties.
function noisecalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to noisecalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addlb.
function addlb_Callback(hObject, eventdata, handles)
% hObject    handle to addlb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename, pathname] = uigetfile('*.mat', 'Additional Library File');
if strcmp(num2str(filename),'0')==0
    set(hObject,'string',[pathname filename]);
else
    set(hObject,'string','-');
end

% --- Executes on button press in process.filedirectory.
function FileDirectory_Callback(hObject, eventdata, handles)
% hObject    handle to process.filedirectory (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.FileDirectory,'string',uigetdir('C:/Jason/Data', 'Select File Directory'))


% --- Executes during object creation, after setting all properties.
function initprocessstatus_CreateFcn(hObject, eventdata, handles)
% hObject    handle to initprocessstatus (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


%no longer using shot noise CM 12/13/2021
%bookmark CM 12/13/2021
%function [outspec,outnoise,outempnoise,cosray] = CRRPolyFit(wavelength,spec,shotnoise,order,c1,c2,mediandeviations,addlb)
function [outspec,outempnoise,cosray] = CRRPolyFit(wavelength,spec,order,c1,c2,mediandeviations,addlb)

% Initialize variables
if iscell(spec) == 0
    spect = spec; clear spec;
    spec{1} = spect; clear spect;
    
    %shotnoiset = shotnoise; clear shotnoise; %CM 12/13/2021
    %shotnoise{1} = shotnoiset; clear shotnoiset; %CM 12/13/2021
end
outspec = cell(size(spec)); outnoise = outspec;

for ijk = 1:size(spec,2)
    for jkl = 1:size(spec{ijk},2)
        rawspect = squeeze(spec{ijk}(:,jkl,:));
        %shott = squeeze(shotnoise{ijk}(:,jkl,:)); %CM 12/13/2021
        [spect] = IanitaJ(wavelength,rawspect,order,c1,c2,1,addlb);
        
        bkgd = rawspect-spect;
        mspec = repmat(median(spect,2),1,size(spect,2));
        medad = repmat(median(abs(spect-mspec),2),1,size(spect,2));
        cutoff = mspec+mediandeviations.*medad;
        
        % Determine indices of cosmic rays
        ind = spect>cutoff;
        ind(1:c1,:) = 0; ind(c2:end,:) = 0;
        
        cosray(ijk,jkl) = sum(ind(:));
        
        % Set cosmic rays to NaN
        spect(ind) = nan; %shott(ind) = nan; %CM 12/13/2021
        
        % sum frames
        outspect = size(spect,2).*nanmean(spect,2);
        outspec{ijk}(:,jkl) = outspect + sum(bkgd,2);
        %outnoise{ijk}(:,jkl) = sqrt(size(spect,2).*nanmean(shott.^2,2));
        %CM 12/13/2021
        
        outempnoise{ijk}(:,jkl) = sqrt(size(spect,2).*nanvar(spect,[],2));
    end
end



function bkgdpolyorder_Callback(hObject, eventdata, handles)
% hObject    handle to bkgdpolyorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bkgdpolyorder as text
%        str2double(get(hObject,'String')) returns contents of bkgdpolyorder as a double



function smoothpolyorder_Callback(hObject, eventdata, handles)
% hObject    handle to smoothpolyorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of smoothpolyorder as text
%        str2double(get(hObject,'String')) returns contents of smoothpolyorder as a double

% --- Executes during object creation, after setting all properties.
function smoothpolyorder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to smoothpolyorder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function [peakind,peakheight] = PeakLocationsJ(xaxis,spec,peaknum,peakstripwindow,edgedist1,fitflag,dwin)
% function [peakind,peakheight] = PeakLocationsJ(xaxis,spec,peaknum,peakstripwindow,edgedist,fitflag)
% Pass empty xaxis to simply use pixel number.  peakstripwindow and
% edgedist are always in PIXELS independent of xaxis.
% Modified by Guanping Feng 03/12/2016
if exist('dwin','var') == 0 
    dwin = 1;
end
if exist('edgedist1','var') == 0 || peakstripwindow>edgedist1
    edgedist1 = peakstripwindow;
end
if exist('fitflag','var') == 0
    fitflag = 1;
end
if isempty(xaxis) == 1
    xaxis = (1:size(spec,1)).';
elseif size(xaxis,1) ~= size(spec,1)
    warning('xaxis is not the same size as the length of spec.  Setting xaxis equal to pixel number...')
    xaxis = 1:size(spec,1);
end
relzero = min(spec);

edgedist2 = (size(spec,1)-edgedist1);
ijk = 0; peakind = zeros(peaknum,1); peakheight = peakind;
while ijk<=peaknum-1
    ijk
    [~,peakindt] = max(spec);
    if peakindt <= edgedist1 || peakindt >= edgedist2
        spec(peakindt) = relzero;
    else
        peakwindow = (peakindt + (-peakstripwindow:peakstripwindow)).';
        if fitflag == 1
            [coeff,a,goodness] = GLFitJ(xaxis(peakwindow),spec(peakwindow),1); %%Find peaks with Guassian fitting
            ijk=ijk+1;
            peakind(ijk) = coeff.x01;
            peakheight(ijk) = coeff.H1;
            [~,tind]=max(spec(peakwindow));
            peakpos=peakwindow(tind);
        else
            [pkh,tind] = max(spec(peakwindow));  %%Find peaks with maxima Seeking
            peakpos=peakwindow(tind);
            if  pkh>spec(peakpos-2)&&pkh>spec(peakpos+2)
                ijk=ijk+1;
                if ijk>peaknum
                    break
                end
                peakheight(ijk)=pkh;
                peakind(ijk) = xaxis(peakwindow(tind));
            end
        end
        %Remove old peak which has been already found
        for iii=peakpos :length(spec)-1
            if  iii+dwin>=length(spec)
                spec(iii:end)=relzero;
                break
            elseif spec(iii)-spec(iii+dwin)>0
                spec(iii)=relzero;
            else
                break
            end
        end
        for iii=peakpos-1:-1:2
            if  iii-dwin<=1
                spec(1:iii)=relzero;
                break
            elseif spec(iii)-spec(iii-dwin)>0
                spec(iii)=relzero;
            else
                break           
            end
        end
    end
end
[peakind,ind] = sort(peakind);
peakheight = peakheight(ind);





function [c,OLSfit,components,residual] = OLSJ(spec,basis)
% function [c,OLSfit,components,residual] = OLSJ(spec,basis)
% Perform an ordinary least squares fit to a set of spectra with a matrix
% composed of the basis spectra for the fit.  The fit coefficients are
% returned in the matrix c.  OLSfit is the fit to each spectrum, and
% components is a cell array with each cell containing a matrix of the
% component spectra that make up the total fit.  All spectra are defined as
% COLUMN vectors.

% Not sure why I can't use basis instead of bT' below, but I get
% strange results when fitting small spectral regions if basis is used
bT = basis';
c = (bT*(bT'))\(bT*spec); % The fit coefficients
% c = basis\spec;
OLSfit = basis*c; % The computed fit
components = cell(1,size(c,2)); % Components that make up the fit
for ijk = 1:size(c,2)
    components{ijk} = basis*diag(c(:,ijk));
end
residual = spec-OLSfit; % Residual

function [specfit,area,areai] = GLEvalJ(x,coeff)
% Evaluate the Gaussian Lorentzian function with coefficients in coeff at
% the locations given in the vector x.  Also returns tbe area under the
% Gaussian Lorentzian curve.

tval = coeffvalues(coeff);
npeaks = round(length(tval)./4);

H = repmat(tval(1:npeaks),length(x),1);
M = repmat(tval(npeaks+1:2*npeaks),length(x),1);
w = repmat(tval(2*npeaks+1:3*npeaks),length(x),1);
x0 = repmat(tval(3*npeaks+1:4*npeaks),length(x),1);
x = repmat(x,1,npeaks);
specfit = (1-M).*H.*exp(-((x-x0)./w).^2.*4.*log(2)) + ...
    M.*H./(4.*((x-x0)./w).^2 + 1);

areas = sprintf(['(1-M(1,%d)).*H(1,%d).*4.*sqrt(pi).*w(1,%d).*log(2) + '...
    'M(1,%d).*H(1,%d).*pi.*w(1,%d)./2'],ones(1,6));
areai = zeros(npeaks,1); % Preallocate
areai(1) = eval(areas);
for ijk = 2:npeaks
    areas = ([areas sprintf([' + (1-M(1,%d)).*H(1,%d).*4.*sqrt(pi).*w(1,%d).*log(2) + '...
        'M(1,%d).*H(1,%d).*pi.*w(1,%d)./2'],ones(1,6).*ijk)]);
    areai(ijk) = eval(sprintf([' + (1-M(1,%d)).*H(1,%d).*4.*sqrt(pi).*w(1,%d).*log(2) + '...
        'M(1,%d).*H(1,%d).*pi.*w(1,%d)./2'],ones(1,6).*ijk));
end

area = eval(areas);

function [coeff,area,goodness,output] = GLFitJ(x,spec,npeaks)
% Performs a mixed Gaussian/Lorentzian Fit of a spectral region.  The
% Gaussian/Lorentzian fit is returned as a funciton of x, the spectral
% region to fit is spec, and the number of peaks used in the fit is npeaks.

% M = mixture (% Lorentzian) between 0 and 1
maxeval = 600E2; % 600 by default
functiontol = 1e-8; % 1e-6 by default
coeftol = 1e-8; % 1e-6 by default

% maxeval = 600; % 600 by default
% functiontol = 1e-6; % 1e-6 by default
% coeftol = 1e-6; % 1e-6 by default

stepsize = floor(length(x)./npeaks);
[he] = max(spec); ce = (x(1:stepsize:length(x))).'; ce = ce(1:npeaks);
wi = x(end)-x(1); Ms = 0.5; he = he./npeaks;

t = ones(1,npeaks);
fitinit = [t.*he t.*Ms t.*wi./2 t.*ce];
lower = [t.*0 t.*0 t.*0 t.*x(1)];
upper = [t.*3.*he.*npeaks t.*1 t.*wi t.*x(end)];

% Mixed Gaussian/Lorentzian Function (as defined by Thermo Scientific)
fn = sprintf(['(1-M%d).*H%d.*exp(-((x-x0%d)./w%d).^2.*4.*log(2)) + '...
    'M%d.*H%d./(4.*((x-x0%d)./w%d).^2 + 1)'],ones(1,8));
for ijk = 2:npeaks
    fn = ([fn ' + ' sprintf(['(1-M%d).*H%d.*exp(-((x-x0%d)./w%d).^2.*4.*log(2)) + '...
        'M%d.*H%d./(4.*((x-x0%d)./w%d).^2 + 1)'],ones(1,8).*ijk)]);
end

fitopt = fitoptions('Method','NonlinearLeastSquares','Startpoint',fitinit,...
    'Lower',lower,'Upper',upper,'MaxFunEvals',maxeval,'TolFun',functiontol,...
    'TolX',coeftol);
f = fittype(fn,'options',fitopt);
% Perform the Non Linear Least Squares Fit
[coeff,goodness,output] = fit(x,spec,f);
[specfit,area] = GLEvalJ(x,coeff);

% Optional Plotting Section
ploton = 0;
if ploton == 1
    figure(99); clf; hold on
    plot(x,spec)
    plot(coeff)
    plot(x,specfit,'k:')
    title(['Rsquared = ' num2str(goodness.rsquare) ...
        '  ::  Area under curve: ' num2str(area)])
    legend('Data','Total Fit','Fit Components')
    pause;
end


function [coeff,area,goodness,output,basis] = GFitFiberBasisJ(x,spec,npeaks,ce,he)
% Performs a Gaussian Lorentzian Fit of a spectral region.  The
% Gaussian fit is returned as a funciton of x, the spectral
% region to fit is spec, and the number of peaks used in the fit is npeaks.

% M = mixture (% Lorentzian) between 0 and 1
maxeval = 10000; % 600 by default
functiontol = 1e-11; % 1e-6 by default
coeftol = 1e-11; % 1e-6 by default

maxeval = 600; % 600 by default
functiontol = 1e-6; % 1e-6 by default
coeftol = 1e-6; % 1e-6 by default


% Temporary
ce = ce.'; he = he.';
t = ones(1,npeaks);
fitinit = [t.*he t.*3.5 t.*ce];
lower = [t.*0 t.*2.5 t.*(ce-1)];
upper = [t.*3.*he t.*6 t.*(ce+1)];

% Note: Never chooses Lorentzian if I do mixed Gaussian Lorentzian
% Gaussian Function
fn = sprintf('H%d.*exp(-((x-x0%d)./w%d).^2.*4.*log(2))',ones(1,3));
for ijk = 2:npeaks
    fn = ([fn ' + ' sprintf('H%d.*exp(-((x-x0%d)./w%d).^2.*4.*log(2))',ones(1,3).*ijk)]);
end

fitopt = fitoptions('Method','NonlinearLeastSquares','Startpoint',fitinit,...
    'Lower',lower,'Upper',upper,'MaxFunEvals',maxeval,'TolFun',functiontol,...
    'TolX',coeftol);
f = fittype(fn,'options',fitopt);
% Perform the Non Linear Least Squares Fit
[coeff,goodness,output] = fit(x,spec,f);
[basis,area,areai] = GEvalJ(x,coeff);
basis = fliplr(basis);

% Optional Plotting Section
ploton = 0;
if ploton == 1
    figure(99); clf; hold on
    plot(x,spec)
    plot(coeff)
    plot(x,basis','k:')
    title(['R^{2} = ' num2str(goodness.rsquare) ...
        '  ::  Area under curve: ' num2str(area)])
    legend('Data','Total Fit','Fit Components')
    pause;
end

function [specfit,area,areai] = GEvalJ(x,coeff)
% Evaluate the Gaussian function with coefficients in coeff at
% the locations given in the vector x.  Also returns tbe area under the
% Gaussian curve.  areai is the area of each individual peak.

tval = coeffvalues(coeff);
npeaks = round(length(tval)./3);

H = repmat(tval(1:npeaks),length(x),1);
w = repmat(tval(1*npeaks+1:2*npeaks),length(x),1);
x0 = repmat(tval(2*npeaks+1:3*npeaks),length(x),1);
x = repmat(x,1,npeaks);
specfit = H.*exp(-((x-x0)./w).^2.*4.*log(2));

areas = sprintf('H(1,%d).*4.*sqrt(pi).*w(1,%d).*log(2)',ones(1,2));
areai = zeros(npeaks,1); % Preallocate
areai(1) = eval(areas);
for ijk = 2:npeaks
    areas = ([areas sprintf(' + H(1,%d).*4.*sqrt(pi).*w(1,%d).*log(2)',ones(1,2).*ijk)]);
    areai(ijk) = eval(sprintf('H(1,%d).*4.*sqrt(pi).*w(1,%d).*log(2)',ones(1,2).*ijk));
end

area = eval(areas);


function [hours, minutes, seconds] = RemainingTimeJ(curriter,totiter)
t = toc;
remainiter = totiter-curriter;
t = t.*remainiter./curriter;
hours = floor(t./3600); t = t-hours.*3600;
minutes = floor(t./60); t = t-minutes.*60;
seconds = round(t);
fprintf([num2str(curriter) ' of ' num2str(totiter) ' - Remaining Time: ' num2str(hours) ' hours ' ...
    num2str(minutes) ' minutes ' num2str(seconds) ...
    ' seconds\n']);


function [spec,components] = IanitaJ(wavelength,spec,order,c1,c2,iter,addlb,convergepercent,peakremovalflag,noisefac,noise)
% Improved automated fluorescence subtraction from Raman spectra
%
% [fitspec,components] = IanitaJ(wavelength,spec,order,c1,c2,iter)
% fitspec: spectrum after fluorescence background subtraction;
% wavelength: wavelength values of measured spectrum; spec: input spectrum;
% order: polynomial order for the fit; c1, c2: fitting region;
% iter: maximum number of iterations
%
% [fitspec,components] = IanitaJ(wavelength,spec,order,c1,c2,iter,addlb)
% addlb: additional spectra to include in the basis set for the fit, use []
% for no additional spectra
%
% [fitspec,components] = IanitaJ(wavelength,spec,order,c1,c2,iter,addlb,convergepercent,peakremovalflag,noisefac,noise)
% convergepercent: stops iterating if percent change in standard deviation
% of residual is less than converge percent; peakremovalflag: set to 1 in
% order to attempt removal of Raman peaks before fitting; noisefac: number
% of standard deviations in noise to permit below zero intensity
% noise: noise at each pixel (same dimension as spec), if noise is not
% passed, IanitaJ makes an estimate of the noise based upon the
% standard deviation of the residual of the fit
%
% Defaults: addlb = []; convergepercent = 5E-4; peakremovalflag = 0; noisefac = 0;
%
% Modeled after 'Automated Autofluorescence Background Subtraction Algorithm
% for Biomedical Raman Spectroscopy' by Zhao et al.
% JRM 11/19/10

% Number of standard deviations in noise to permit below zero intensity
if exist('noisefac','var') == 0
    noisefac = 0; % noisefac = 0 performs a conventional ANITA fit
end
if exist('peakremovalflag','var')==0
    peakremovalflag = 0; % peakremovalflag = 1 attempts to remove Raman peaks before fitting
end
% Convergence Condition
if exist('convergepercent','var') == 0
    convergepercent = 5E-4;
end
if exist('noise','var') == 1
    noiseflag = 1;
else
    noiseflag = 0;
end

% Create Polynomials to Use in Basis
px = (wavelength-mean(wavelength))./std(wavelength);
% Initialize polynomial matrix
polyb = zeros(size(spec,1),order+1);
% Values range from 1 to 2^order
for j=0:1:order
    polyb(:,j+1)=px.^j;
end

% Add additional Basis Vectors
if exist('addlb','var') == 1
    basis = [polyb addlb];
else
    basis = polyb;
end

% For each spectrum
for jkl = 1:size(spec,2)
    putin = spec(c1:c2,jkl); reducedbasis = basis(c1:c2,:);
    
    % Perform the Ordinary Least Squares Fit Between c1 and c2
    [~,putout,~,residual] = OLSJ(putin,reducedbasis);
    DEV = std(residual);
    
    % Raman peak removal
    if noiseflag == 1
        noisecorrection = noise(c1:c2,jkl);
        if peakremovalflag == 1
            noremoval = (putin<(putout+3.*noisecorrection));
            noisecorrection = noisecorrection(noremoval);
        end
    elseif peakremovalflag == 1
        noremoval = (putin<(putout+3.*DEV));
    end
    if peakremovalflag == 1
        reducedbasis = reducedbasis(noremoval,:);
        putin = putin(noremoval);
        if isreal(putin) == 0
            pause;
        end
    end
    
    % 1st iteration
    % Perform the Ordinary Least Squares Fit Between c1 and c2
    [c,~,components,residual] = OLSJ(putin,reducedbasis);
    DEVlast = std(residual);
    
    if noiseflag ~= 1 % Estimate noise as DEV
        noisecorrection = DEVlast;
    end
    
    % Reasign fit values that would result in negative spectral
    % intensity
    residualcorrection = residual-noisefac.*noisecorrection;
    cutoff = (residualcorrection+abs(residualcorrection))./2;
    
    putin = putin - cutoff;
    
    for ijk = 2:iter % For each iteration
        % Perform the Ordinary Least Squares Fit Between c1 and c2
        [c,~,components,residual] = OLSJ(putin,reducedbasis);
        
        % Calculate the standard deviation of the residual
        DEV = std(residual);
        
        differ = abs((DEV-DEVlast)./DEV);
        % If the difference is small break from the loop
        if differ < convergepercent
            break
        end
        DEVlast = DEV;
        
        if noiseflag ~= 1 % Estimate noise as DEV
            noisecorrection = DEV;
        end
        
        % Reasign fit values that would result in negative spectral
        % intensity
        residualcorrection = residual-noisefac.*noisecorrection;
        cutoff = (residualcorrection+abs(residualcorrection))./2;
        
        putin = putin - cutoff;
    end % of for each iteration
    
    % Subtract estimated fluorescence lineshape
    spec(:,jkl) = spec(:,jkl) - basis*c;
end % of for each spectrum


function Smooth = SvgSmoothJ(M,n,Spec)
% Author: Jeff Rice
% Date: 7-17-02
% Smooth a given spectrum using a Sovitsky-Golay filter
% This function accepts an input spectrum and applies a Sovitsky-Golay filter to it,
% returning the smoothed spectrum
% "M" is the order of the Savitsky-Golay filter
% "n" is the number of data points fit on each side of the point of interest
% i.e. the total number of points fit is 2n+1


Spec = Spec'; %JRM

A1 = [-n:n];                % create a row vector with values -n,-n+1,....n
%     A1 = A1./std(A1); normalize as well??
A = ones(M+1,2*n+1);        % actually it's A transpose.

for j = 1:M,
    A(j+1,:) = A(j,:).*A1;  % fill the second row of A with A1*A1 or A1^2, etc.
end

P = A*A';                   % If A is A(but not its transpose), here it should be P = A'*A;
Inv = inv(P);

% c= zeros(1,n);
C = zeros(1,2*n+1);         % form C, a row vector to store coefficients

for I = 1:n;                % form I, a vector of n^m where m = 0,1,...M
    K = ones(M+1,1);
    
    for J = 1:M,
        K(J+1) = K(J)*I;
    end
    
    %     c(I) = Inv(1,:)*K;
    c = Inv(1,:)*K;
    C(n+1+I) = c;
    C(n+1-I) = c;
end

C(n+1) = Inv(1,1);

Spec1 = [zeros(1,n)+Spec(1) Spec zeros(1,n)+Spec(end)];
Smooth = conv(Spec1,C);
Smooth = Smooth(2*n+1:(end-2*n));
Smooth = Smooth'; %JRM


function [cs,index] = sort_nat(c)
%sort_nat: Natural order sort of cell array of strings.
% usage:  [S,INDEX] = sort_nat(C)
%
% where,
%    C is a cell array (vector) of strings to be sorted.
%    S is C, sorted in natural order.
%    INDEX is the sort order such that S = C(INDEX);
%
% Natural order sorting sorts strings containing digits in a way such that
% the numerical value of the digits is taken into account.  It is
% especially useful for sorting file names containing index numbers with
% different numbers of digits.  Often, people will use leading zeros to get
% the right sort order, but with this function you don't have to do that.
% For example, if C = {'file1.txt','file2.txt','file10.txt'}, a normal sort
% will give you
%
%       {'file1.txt'  'file10.txt'  'file2.txt'}
%
% whereas, sort_nat will give you
%
%       {'file1.txt'  'file2.txt'  'file10.txt'}
%
% See also: sort

% Version: 1.2, 5 November 2008
% Author:  Douglas M. Schwarz
% Email:   dmschwarz=ieee*org, dmschwarz=urgrad*rochester*edu
% Real_email = regexprep(Email,{'=','*'},{'@','.'})


% Replace runs of digits with '0'.
c2 = regexprep(c,'\d+','0');

% Compute char version of c2 and locations of zeros.
s1 = char(c2);
z = s1 == '0';

% Extract the runs of digits and their start and end indices.
[digruns,first,last] = regexp(c,'\d+','match','start','end');

% Create matrix of numerical values of runs of digits and a matrix of the
% number of digits in each run.
num_str = length(c);
max_len = size(s1,2);
num_val = NaN(num_str,max_len);
num_dig = NaN(num_str,max_len);
for i = 1:num_str
    num_val(i,z(i,:)) = sscanf(sprintf('%s ',digruns{i}{:}),'%f');
    num_dig(i,z(i,:)) = last{i} - first{i} + 1;
end

% Find columns that have at least one non-NaN.  Make sure activecols is a
% 1-by-n vector even if n = 0.
activecols = reshape(find(~all(isnan(num_val))),1,[]);
n = length(activecols);

% Compute which columns in the composite matrix get the numbers.
numcols = activecols + (1:2:2*n);

% Compute which columns in the composite matrix get the number of digits.
ndigcols = numcols + 1;

% Compute which columns in the composite matrix get chars.
charcols = true(1,max_len + 2*n);
charcols(numcols) = false;
charcols(ndigcols) = false;

% Create and fill composite matrix, comp.
comp = zeros(num_str,max_len + 2*n);
comp(:,charcols) = double(s1);
comp(:,numcols) = num_val(:,activecols);
comp(:,ndigcols) = num_dig(:,activecols);

% Sort rows of composite matrix and use index to sort c.
[unused,index] = sortrows(comp);
cs = c(index);


function h = plotR(x,y,errorpatch,varargin)
% function plotR(x,y,errorpatch,varargin)
% Plots a Raman spectrum with common labels

% Defaults
fontsizelabel = 8;
fontsizeaxis = 8;
fontweight = 'bold';

h = plot(x,y,varargin{:},'DisplayName','y');
if exist('errorpatch','var') == 1
    if isempty(errorpatch) == 0
        if size(errorpatch,3)>1
            errorpatchL = squeeze(errorpatch(:,:,1));
            errorpatchU = squeeze(errorpatch(:,:,2));
        else
            errorpatchL = errorpatch;
            errorpatchU = errorpatch;
        end
        
        xp = [x; flipud(x)];
        for ijk = 1:size(h,1)
            linecol = get(h(ijk),'color');
            
            delete(h(ijk));
            yp = [y(:,ijk) + errorpatchU(:,ijk); flipud(y(:,ijk) - errorpatchL(:,ijk))];
            
            linecol = linecol+0.75; linecol(linecol>1) = 1;
            patch(xp,yp,linecol,'EdgeColor','none'); hold on;
            
            %             patch(xp,yp,linecol,'EdgeColor','none'); hold on
        end
        h = plot(x,y,varargin{:},'DisplayName','y');
    end
    
end
xlabel('Raman shift / cm^{-1}','fontsize',fontsizelabel,'fontweight',fontweight)
ylabel('Raman intensity / a.u.','fontsize',fontsizelabel,'fontweight',fontweight)
set(gca,'Ytick',[],'fontsize',fontsizeaxis,'layer','top','linewidth',0.75)
axis tight; box off;
yt = ylim; yrange = range(yt);
yt(1) = yt(1)-yrange.*0.05; yt(2) = yt(2)+yrange.*0.05; ylim(yt);
