function varargout = Settings_GUI(varargin)
% AUTHOR:   Abel Szkalisity
% DATE:     Oct 20, 2019
% NAME:     Settings_GUI
%
%      General settings GUI. Designed for ask parameters from users. The
%      actual parameters are described in the inputs.
%      This GUI is NOT singleton, but check in GUIDE.
%      Use the arguments as specified here.
%
% INPUTS:
%    paramarray     A cellarray with structs that describe the
%                   parameters. See doc of generateUIControls.
%
%  Optional NAME-VALUE pairs:
%    infoText       A string text to be displayed as an info text at the
%                   top of the GUI.
%    infoAlignment  A string specifying how to align the info text
%                   horizontally. Must be one of {'center', 'left',
%                   'right'} (default: left)
%    title          The title of the figure to be displayed on the main bar
%                   (NOTE: some Operating systems might not display this)
%    checkFunction  Function handle to the evaluation function which can be
%                   used to check if the parameters are fulfilling some
%                   requirements. Please see the checkHandle documentation
%                   of the fetchUIControlValues.m file.
%   checkVisibility See the documentation in generateUIControls.
%   checkInterimValues boolean. The dynamic visualization feature provided
%                   by the checkVisibility function handle above may be
%                   easier to implement if the visibility function does not
%                   have to check whether all parameters provided by the
%                   user are faithful to the description in paramarray
%                   (e.g. numbers are really numbers or just some arbitrary
%                   string). Generally parameters are checked to be correct
%                   when the GUI returns with the answers, but not whilst
%                   the user is still in the process of filling.
%                   Consequently by default, the checkVisibility function's
%                   paramcell input is NOT checked by the provided
%                   checkFunction (or the default checks if it is not
%                   provided). This functionality is encoded in the default
%                   false value of this boolean variable. Set it to true if
%                   you want to check also the parameters before sending it
%                   to the checkVisibility function. In this case you do
%                   not need to check the integrity of the returned values
%                   within the checkVisibility function itself. NOTE: if
%                   you set this to true, you must pay attention to proper
%                   default variables, as in that case the check functions
%                   are also called upon the first loading. In case of
%                   failure every input field will be shown.
%          
%
% OUTPUT:
%    paramcell      A cellarray with equal length to paramarray. It
%                   contains the user specified entries in the format
%                   matching the corresponding paramarray entry and as
%                   written in fetchUIControlValues.
%
%    NOTE: in case the user press the close button this function raises an
%    error with the following ID: 'Settings_GUI:noParameterProvided'
%
%   See also generateUIControls, fetchUIControlValues, fillOutUIControls
%
% COPYRIGHT
% Settings Template Toolbox. All Rights Reversed. 
% Copyright (C) 2019 Abel Szkalisity
% BIOMAG, Synthetic and System Biology Unit, Institute of Biochemsitry,
% Biological Research Center, Szeged, Hungary
% Ikonen group Department of Anatomy, Faculty of Medicine, University of
% Helsinki, Helsinki, Finland.
%
% Last Modified by GUIDE v2.5 20-Oct-2019 08:49:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Settings_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Settings_GUI_OutputFcn, ...
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


%Creating extra panel for the scrollbar option
% --- Executes just before Settings_GUI is made visible.
function Settings_GUI_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Settings_GUI (see VARARGIN)

% Choose default command line output for Settings_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

uData = get(handles.figure1,'UserData');

% Parsing inputs
p = inputParser;
addRequired(p,'paramarray');
addParameter(p,'infoText','Specify parameters:');
addParameter(p,'checkFunction',[]);

validAlignments = {'center','left','right'};
checkAlignment = @(x) any(validatestring(x,validAlignments));
addParameter(p,'infoAlignment','left',checkAlignment);
addParameter(p,'title','Settings Dialog',@ischar);
addParameter(p,'checkVisibility',@(a,b)(defaultCheckVisibility(a,b,true)));
addParameter(p,'checkInterimValues',false);
parse(p,varargin{:});

paramarray = p.Results.paramarray;
infoText = p.Results.infoText;

if ~isempty(p.Results.checkFunction)
    uData.checkFunction = p.Results.checkFunction;
end
uData.checkVisibility = p.Results.checkVisibility;

% Set parental relation for the scrollbar
set(handles.uiPanel_InnerSettings,'Parent',handles.uiPanel_Settings);

%Set up GUI default sizes
uData = setUpGUIsize(uData);

set(handles.text_Info,'String',infoText);
if p.Results.checkInterimValues
[ uData.paramTextHandles, uData.paramEditHandles] = generateUIControls(...
    paramarray, handles.uiPanel_InnerSettings, uData.sizeParameters.itemHeight, uData.sizeParameters.margin,...
    'resizeParent',true,...
    'divisionRatio',uData.sizeParameters.textRatio,...
    'checkVisibility',p.Results.checkVisibility,...
    'checkHandle',p.Results.checkFunction);
else
    [ uData.paramTextHandles, uData.paramEditHandles] = generateUIControls(...
    paramarray, handles.uiPanel_InnerSettings, uData.sizeParameters.itemHeight, uData.sizeParameters.margin,...
    'resizeParent',true,...
    'divisionRatio',uData.sizeParameters.textRatio,...
    'checkVisibility',p.Results.checkVisibility);
end
uData.paramarray = paramarray;

%Calculate default opening height and some sizeParameters:
stringLineWidth = uData.sizeParameters.stringLineWidth;
    
%Calculate top size
infoString = get(handles.text_Info,'String');
nofLines = size(infoString,1);
uData.sizeParameters.infoHeight = nofLines*stringLineWidth;

defaultOpenHeight = handles.uiPanel_InnerSettings.Position(4) + uData.sizeParameters.infoHeight + 2*uData.sizeParameters.innerMargin + uData.sizeParameters.buttonHeight + 2*uData.sizeParameters.sideMargin;
handles.figure1.Position(4) = min(defaultOpenHeight,uData.sizeParameters.maxDefaultFigureHeight)+20; %+20 magic const to avoid the scrollbar on first start

% Set alignment accordingly
set(handles.text_Info,'HorizontalAlignment',p.Results.infoAlignment);
set(handles.figure1,'Name',p.Results.title);

% START of setting scrollbar
hPanel = handles.uiPanel_InnerSettings;
% Get the panel's underlying JPanel object reference
jPanel = hPanel.JavaFrame.getGUIDEView.getParent;
 
% Embed the JPanel within a new JScrollPanel object
jScrollPanel = javaObjectEDT(javax.swing.JScrollPane(jPanel));
 
% Remove the JScrollPane border-line
jScrollPanel.setBorder([]);
 
% Place the JScrollPanel in same GUI location as the original panel
pixelpos = getpixelposition(hPanel);
hParent = hPanel.Parent;
[hjScrollPanel, hScrollPanel] = javacomponent(jScrollPanel, pixelpos, hParent);
hScrollPanel.Units = 'norm';
hScrollPanel.OuterPosition = [0 0 1 1];

% Ensure that the scroll-panel and contained panel have linked visibility
hLink = linkprop([hPanel,hScrollPanel],'Visible');
setappdata(hPanel,'ScrollPanelVisibilityLink',hLink);
%END of custom scroll panel

set(handles.figure1,'UserData',uData);

set(handles.figure1,'Visible','on');
pause(0.05); %this is needed so that the setValue works properly
hjScrollPanel.getVerticalScrollBar().setValue(0);
hjScrollPanel.getVerticalScrollBar().setUnitIncrement(5)

% UIWAIT makes Settings_GUI wait for user response (see UIRESUME)
uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Settings_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isfield(handles,'figure1')
    uData = get(handles.figure1,'UserData');
    varargout{1} = uData.paramcell;
    figure1_CloseRequestFcn(hObject, eventdata, handles);
else
    error('Settings_GUI:noParameterProvided','The required parameters were not closed.');
    %varargout{1} = [];
end


% --- Executes on button press in Button_OK.
function Button_OK_Callback(~, ~, handles) %#ok<DEFNU> Called from GUI
% hObject    handle to Button_OK (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uData = get(handles.figure1,'UserData');
try
    if isfield(uData,'checkFunction')
        [uData.paramcell] = fetchUIControlValues( uData.paramarray,uData.paramEditHandles,uData.checkFunction);
    else
        [uData.paramcell] = fetchUIControlValues( uData.paramarray,uData.paramEditHandles);
    end
catch ME
    %if something is wrong with the parameters then we return
    if strcmp(ME.identifier,'Settings_GUI:errorConstraintFault')
        return;
    else
        rethrow(ME);
    end
end
set(handles.figure1,'UserData',uData);
uiresume(handles.figure1);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, ~, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(handles.figure1);
delete(hObject);


% --- Executes when figure1 is resized.
function figure1_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Parameters
uData = get(handles.figure1,'UserData');
sideMargin = uData.sizeParameters.sideMargin;
innerMargin = uData.sizeParameters.innerMargin;
buttonWidth = uData.sizeParameters.buttonWidth;
buttonHeight = uData.sizeParameters.buttonHeight;
infoHeight = uData.sizeParameters.infoHeight;
scrollBarSize = uData.sizeParameters.scrollBarSide;

figurePos = get(hObject,'Position');
set(handles.text_Info,'Position',[sideMargin,figurePos(4)-sideMargin-infoHeight,figurePos(3)-2*sideMargin,infoHeight]);
set(handles.uiPanel_Settings,'Position',[sideMargin,innerMargin+sideMargin+buttonHeight,figurePos(3)-2*sideMargin,figurePos(4)-2*innerMargin-2*sideMargin-infoHeight-buttonHeight]);
requiredHeight = handles.uiPanel_InnerSettings.Position(4);
set(handles.uiPanel_InnerSettings,'Position',[scrollBarSize/2,0,figurePos(3)-2*sideMargin-scrollBarSize ,requiredHeight]);
set(handles.Button_OK,'Position',[(figurePos(3)-buttonWidth)/2,sideMargin,buttonWidth,buttonHeight]);

placeUIControls( ...
    uData.paramarray, uData.paramTextHandles, uData.paramEditHandles,...
    handles.uiPanel_InnerSettings.Position, uData.sizeParameters.itemHeight, uData.sizeParameters.margin,...
    'divisionRatio',uData.sizeParameters.textRatio,...
    'checkVisibility',uData.checkVisibility);
    
function uData = setUpGUIsize(uData)

uData.sizeParameters.stringLineWidth = 15;
uData.sizeParameters.sideMargin = 20;
uData.sizeParameters.innerMargin = 20;
uData.sizeParameters.buttonWidth = 100;
uData.sizeParameters.buttonHeight = 25;
uData.sizeParameters.maxDefaultFigureHeight = 650;
uData.sizeParameters.itemHeight = 20;
uData.sizeParameters.margin = [10 10];
uData.sizeParameters.scrollBarSide = 25;
uData.sizeParameters.textRatio = 0.5;

