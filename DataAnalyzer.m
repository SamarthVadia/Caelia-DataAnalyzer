%% Content
%{
Functions in UI in order they are defined
icon_update() - icons for various buttons are defined
updateimgidlist() - update imageIDs for viewing list (continuously updated) (change number of images in list here)
updatefitlist() - update list of fitting functions
loadPCAbasis() - loads 'pwoaSave.mat' to create the starting basis for PCA
makeEigenBasis(imageList) - makes an eigenbasis of the images in the argument
updateanalysisdblist() - update imageIDs for analysis list 
updatedblist() - update viewing list with help of imageIDs from updateimgidlist() 
dblist_click(source, eventdata) - load image and info of selected imageID from viewing list(dblist)
dblist_keypress(source, eventdata) - executes del_click on pressing delete key in dblist
updatecurrentimginfo(currentimgid) - update all info boxes for selected imageID in viewing list
framelist_click(source, eventdata) - load image on changing frame (Final, PWA, PWOA, DF)
pca_click(~, ~) - load image on turning on/off PCA mode
hdwmode_change(source, eventdata) - change mode (Normal vs Kinetics)
colormapname_click (~, ~) - change colormap (Color vs B/W)
getImgMode() - gets the image mode requested. 1 = normal, 2 = kinetics with two frames, 3 = kinetics all on one frame
showimg(filenum) - plots image for selected imageID (called by dblist_click and other functions)
doPCA(imgIn) - performs PCA on the selected PWA image using the global eigenBasis
addToBasis(newImg) - adds a new col*row image to the PCA basis and updates the eigenbasis
getImage(imgid,imgmode,framenum) - gets the requested image from the database
data_evaluation(raw_data, imgmode, framenum) - in separate *.m file, determines data from various parameters (for showimg, update, fit functions)
get_roi(imgid,imgmode,framenum) - gets the ROI as set by the cursors
load_click(source, eventdata) - load particular imageID into viewing list
save_click(source, eventdata) - save selected imageID in database permanently
del_click(source, eventdata) - delete selected imageID from database
rename_imgname(source, eventdata) - rename selected imageID in database
add2ana_click(source, eventdata) - add selected imageID(s) in analysis list
analysisdblist_click(source, eventdata) - show imageID of selected image in analysis list
analist_keypress(source, eventdata) - execute del2ana_click in analysis list on pressing delete key
load2ana_click(source, eventdata) - load sequenceID in analysis list
del2ana_click(source, eventdata) - deleted selected imageIDs from analysis list (not from database)
clear2ana_click(source, eventdata) - clear all imageIDs from analysis list (not from database)
zoom_on(source,eventdata) - turn zoom on
zoom_reset(source,~) - zoom reset to 100%
pan_on(source,eventdata) - turn pan on
rotate_on(source,eventdata) - rotate image by given angle
enter_rotangle(source,~) - enter rotation angle
rotinc_click(~,~)/rotdec_click(~,~) - callback function to increment/decrement rotangle
curs_on(source, eventdata) - turn cursor on
curs_update(), xl1_change(s,e), xl2_change(s,e), yl1_change(s,e), yl2_change(s,e) - update cursor position (xcurs, ycurs)
update_but() - updates parameters (n count, width, etc) from image in roi
fit_click(source, eventdata) - fits analysis list to selected fit function
singlefit_click(~, ~) - fits single selected image from analysis list to selected fit function
xvardet(anaimgidlist) - determine x variable for fitting (from file or database)
cftool_click(source, eventdata) - opens cftool with data from results
enter_fitxvar(source, eventdata) - enter position in file name for x variable
xvarinc_click(~,~)/xvardec_click(~,~) - callback function to increment/decrement xvar
yvarinc_click(~,~)/yvardec_click(~,~) - callback function to increment/decrement yvar
updatexvardropmenu() - list of variables from database for x variable
updatefittypelist() - set fit type list (Gaussian, linear, etc) 
getfittype() - determine fit type for fit_click
extractimg_click(~, ~) - to extract specified ROI of image shown in analyzer to tiff, opens prompt to specify destination
deltemp_click(source, eventdata) - delete all images which are stored as 'temp'
closefig_click(source, eventdata) - close all figures except GUI
showanalysisimgid_click(source, eventdata) - show selected image from analysisdblist in img axes
clearlastfit_click(source, eventdata) - clears data saved from last fit (required when fit after changing cursor or angle)
nextimgname_enter(source, eventdata) - To update string of nextimgname edit box on pressing enter
imgname_enter(~, eventdata)- To update string of imgname edit box on pressing enter
globalShortcuts(source, eventdata)- Defines global shortcuts including Update and Fit buttons
onClose(~, ~) - saves PCA data on closing the data analysis instance/figure%}
%}
%%
function DataAnalyzer
% This is version 1.0

%  Create UI 
f = figure('Name', 'Data Analysis Software','Visible','on','Position',[50,50,1600,950],'WindowKeyPressFcn',@globalShortcuts);

[zoom_icon,pan_icon,curs_icon,rotate_icon]=icon_update();

%Database Connection
conn = database('becivdatabase', 'root', 'w0lfg4ng', 'Server', '18.62.27.10', 'Vendor', 'MySQL'); 

% "Global" Variables
col=1024;
row=1024;
xcurs=[1 col];
ycurs=[1 row];
xvar='1';
numImages = 500; % the number of images viewable in the data analyzer
pcaBasisSize = 50;
pwoa = zeros(pcaBasisSize,col*row); % Reshaped probe without atom images to form a basis for pca.
eigenBasis = 0; % Eigenbasis found with PCA.
oldestImgIndex = 1; % Points to the oldest image in pwoaList, which should be replaced first.
imgidlist=[];
anaimgidlist=[];
updateimgidlist();
selecimgidlist=[];
currentimgid=num2cell(max(cell2mat(imgidlist)));
rotangle='0';

%% Define components of UI
img = axes('Units','pixels','Position',[50,330,750,600]);  %Main image from database

% Quick update for data in img
updatebut = uicontrol('Style','pushbutton','String','Update [F5]','Position',[820,550,70,25], 'Callback', @update_but); %Fit files selected from dblist into fitplt
quickres = uicontrol('Style','edit','String','Quick Results','min', 0, 'max', 100, 'Position',[820,380,150,150]); %Quick results for img
autoupbox = uicontrol('Style','checkbox','String','Auto Update','Position',[895,552,120,20], 'Value', 0); %Update the quick results after changing the image.

% Mode of image in img
framelist = uicontrol('Style','listbox', 'min' , 0, 'max' , 1, 'Position', [820, 690, 150,100], 'String', {'Absorption Image','Probe with Atoms','Probe without Atoms','Dark Field','Dark Field (Dual DF)'}, 'Callback', @framelist_click); %List of frames
imgModeText = uicontrol('Style','text','String','Imaging Mode:','Position',[820, 890, 70, 25]);
imgMode = uicontrol('Style','popupmenu','String',{'Normal','Kinetics 2','Kinetics 3'},'Position',[820, 880, 80, 20],'Callback',@hdwmode_change);

% PCA
pcacbox = uicontrol('Style','checkbox','String','Apply PCA','Position',[820,800,90,20],'Value',0,'Callback',@pca_click);
pcashowbox = uicontrol('Style','checkbox','String','Show PCA','Position',[820,825,90,20],'Value',0,'Callback',@pca_click);
pcaTest = uicontrol('Style','checkbox','String','PCA Test','Position',[820,850,90,20],'Value',0,'Callback',@pca_click);

% Tool box for img
zon = uicontrol('Style','togglebutton','CData', zoom_icon,'Position',[855,650,25,25], 'Callback', @zoom_on); %Zoom On for Img
zreset = uicontrol('Style','togglebutton','String', '100%','Position',[890,650,35,25], 'Callback', @zoom_reset); %Zoom On for Img
curs = uicontrol('Style','togglebutton','CData',curs_icon,'Position',[820,650,25,25], 'Callback', @curs_on); %Put cursors in Img axes
panbtn = uicontrol('Style','togglebutton','CData', pan_icon, 'Position',[935,650,25,25], 'Callback', @pan_on); %Turn on pan for Img axes
rotatebtn = uicontrol('Style','togglebutton','CData', rotate_icon, 'Position',[820,610,25,25],'Callback', @rotate_on); %Rotate for Img axes
rotangleinput = uicontrol('Style','edit','String',rotangle,'Position',[855,610,25,25],'Callback', @enter_rotangle); %Angle to rotate image by
colormapname = uicontrol('Style','popupmenu','String',{'Color';'B&W'},'Position',[910,610,50,25], 'Callback', @colormapname_click); %Color map for img axes

rotincrement = uicontrol('Style','pushbutton','String','+','Position',[880,623,25,12],'Callback',@rotinc_click);
rotdecrement = uicontrol('Style','pushbutton','String','-','Position',[880,610,25,12],'Callback',@rotdec_click);

dblist = uicontrol('Style','listbox', 'min' , 0, 'max' , 400, 'Position', [50, 20, 350,280], 'Value', [], 'Callback', @dblist_click, 'KeyPressFcn', @dblist_keypress); %List from database

%Buttons for dblist
savebut = uicontrol('Style','pushbutton','String','Save','Position',[440,145,90,25], 'Callback', @save_click); %Save data from dblist to permanenet database
loadbut = uicontrol('Style','pushbutton','String','Load','Position',[440,200,90,25], 'Callback', @load_click); %Load data manually (from permanenet database)
delbut = uicontrol('Style','pushbutton','String','Delete','Position',[440,90,90,25], 'Callback', @del_click); %Delete data from dblist
add2anabut = uicontrol('Style','pushbutton','String','Add2Analysis(->)','Position',[440,265,90,25], 'Callback', @add2ana_click); %Add files selected in dblist to analysisdblist
addnext = uicontrol('Style','checkbox','String','Add Next','Position',[440,240,70,20], 'Value', 0); %Add the next shot to analysisdblist
nextimgname_text = uicontrol('Style','text','String','Next Image Name: [Ins]','Position',[440,45,120,15]);
nextimgname = uicontrol('Style','edit','Position',[440,20,160,25],'KeyPressFcn', @nextimgname_enter); %Next shot name

%Status boxes for dblist
imgid_text = uicontrol('Style','text','String','Image ID:','Position',[590,275,70,25]);
imgidbox = uicontrol('Style','edit','String','Image ID','Position',[650,280,170,25]); %Display imageID of selected image in dblist
imgname_text = uicontrol('Style','text','String','Image Name:','Position',[580,230,70,25]);
imgname = uicontrol('Style','edit','String','Image Name','Position',[650,235,170,25], 'Callback', @rename_imgname, 'KeyPressFcn', @imgname_enter); %Rename selected image in dblist
timestmp_text = uicontrol('Style','text','String','Timestamp:','Position',[580,185,70,25]);
timestmpstatus = uicontrol('Style','edit','String','Timestamp','Position',[650,190,170,25]); %Displays timestamp of selected image in dblist
savestatus_text = uicontrol('Style','text','String','Saved:','Position',[590,140,37,25]);
savestatus = uicontrol('Style','edit','String','Save Status','Position',[630,145,50,25]); %Displays save status of selected image in dblist
seqidstatus_text = uicontrol('Style','text','String','Sequence ID:','Position',[695,140,70,25]);
seqidstatus = uicontrol('Style','edit','String','ID','Position',[770,145,50,25]); %Displays seq ID of selected image in dblist
seqnamestatus_text = uicontrol('Style','text','String','Sequence Name:','Position',[550,95,90,25]);
seqnamestatus = uicontrol('Style','edit','String','Sequence Name','Position',[650,100,170,25]); %Displays seq name of selected image in dblist
shownamecheck = uicontrol('Style','checkbox','String','Show in list','Position',[690,80,120,20], 'Value', 0); %Show the sequence name in dblist
seqdescstatus_text = uicontrol('Style','text','String','Sequence Description:','Position',[610,40,70,35]);
seqdescstatus = uicontrol('Style','edit','min',0,'max',10,'String','Sequence Description','Position',[680,20,140,60]); %Displays seq description of selected image in dblist

fitplt = axes('Units','pixels','Position',[1050,480,500,450]);  %Fitting function plot
%singplt = axes('Units','pixels','Position',[1280,680,300,250]);  %Single data plot

% Fit button and result
fitbut = uicontrol('Style','pushbutton','String','Fit [F9]','Position',[1455,55,80,30], 'Callback', @fit_click); %Fit files selected from analysisdblist into fitplt
singlefitbut = uicontrol('Style','pushbutton','String','Fit selected image','Position',[1450,20,100,25], 'Callback', @singlefit_click); %Fit single file selected from analysisdblist into fitplt
cftoolbut = uicontrol('Style','pushbutton','String','CFTool','Position',[1240,260,80,25], 'Callback', @cftool_click); %Opens cftool with x & y data from results
fitplotcheck = uicontrol('Style','checkbox','String','Plot fit for each shot','Position',[1320,60,120,20], 'Value', 0); %To check if each shot gives plot for fit
fitres = uicontrol('Style','edit','String','Results','min', 0, 'max', 100,'Position',[1370,260,180,170]); %Results for fit function
analysisimgidbox = uicontrol('Style','edit','String','Image ID','Position',[1280,20,70,25]); %Display imageID of selected image in analysisdblist
showanalysisimgidbut = uicontrol('Style','pushbutton','String','Show','Position',[1360,20,70,25], 'Callback', @showanalysisimgid_click); %Show selected image in analysisdblist
fitoutputnum_text = uicontrol('Style','text','String','Output # :','Position',[1210,295,60,25]);
fitoutputnum = uicontrol('Style','edit','String','1','Position',[1270,300,50,25]); %Number of output variable for fitting

yvarincrement = uicontrol('Style','pushbutton','String','+','Position',[1320,313,25,12],'Callback',@yvarinc_click);
yvardecrement = uicontrol('Style','pushbutton','String','-','Position',[1320,300,25,12],'Callback',@yvardec_click);

% X Variable for fit
xvariable_text = uicontrol('Style','text','String','X Variable:','Position',[1290,220,70,25]);
xvarmode = uibuttongroup('Title','X Var Type','Units', 'pixels', 'Position',[1300, 100, 140,120]); %'SelectionChangedFcn', @hdwmode_change); 
xvarnumbermode = uicontrol('Parent', xvarmode, 'Style','radiobutton','String','None','Units', 'pixels','Position',[10,80,70,25], 'Tag', '1'); %X Variable to be real numbers
xvarnamemode = uicontrol('Parent', xvarmode, 'Style','radiobutton','String','Name','Units', 'pixels','Position',[10,45,70,25], 'Tag', '2'); %Choose x variable for fitplt from name string
fitxvar = uicontrol('Style','edit','String',xvar,'Position',[1450,145,50,25], 'Callback', @enter_fitxvar); %X Variable for fitting function
xvardbmode = uicontrol('Parent', xvarmode, 'Style','radiobutton','String','DB Variable','Units', 'pixels','Position',[10,10,120,25], 'Tag', '3'); %Choose x variable for fitplt from db tables
xvardropmenu = uicontrol('Style','popupmenu','Position',[1450,110,100,25]); %X Variable from database

xvarincrement = uicontrol('Style','pushbutton','String','+','Position',[1500,158,25,12],'Callback',@xvarinc_click);
xvardecrement = uicontrol('Style','pushbutton','String','-','Position',[1500,145,25,12],'Callback',@xvardec_click);

% Y variable for fit
yvariable_text = uicontrol('Style','text','String','Y Variable:','Position',[1050,230,70,25]);
analysisdblist = uicontrol('Style','listbox', 'min' , 0, 'max' , 100, 'Position', [1050, 20, 200,180], 'Callback', @analysisdblist_click, 'KeyPressFcn', @analist_keypress); %List of images to be analyzed
load2anabut = uicontrol('Style','pushbutton','String','Load','Position',[1050,210,50,25], 'Callback', @load2ana_click); %Load images to analysisdb
del2anabut = uicontrol('Style','pushbutton','String','Delete','Position',[1120,210,50,25], 'Callback', @del2ana_click); %Delete images from analysisdb
clear2anabut = uicontrol('Style','pushbutton','String','Clear All','Position',[1190,210,50,25], 'Callback', @clear2ana_click); %Clear all images from analysisdb

% Fit type for fit
fitlist_text = uicontrol('Style','text','String','Measure:','Position',[1040,420,70,15]);
fitlist = uicontrol('Style','listbox', 'min' , 0, 'max' , 1, 'Position', [1050, 280, 140,140]); %List of fitting functions
fittype_text = uicontrol('Style','text','String','Fit with:','Position',[1200,420,70,15]);
fittypelist = uicontrol('Style','listbox', 'min' , 0, 'max' , 1, 'Position', [1210, 335, 130,85]); %List of fit types (gaussian, lorenzian, etc)

% Additional buttons
nextname_text = uicontrol('Style','text','String','Specify next name:','Position',[870,185,100,20]); 
nextnameOption = uicontrol('Style','popupmenu','String',{'None','Text','Text + Variable'},'Position',[870, 165, 120, 20]); %To check how name is specified for new image
extractimgbut = uicontrol('Style','pushbutton','String','Extract roi to tiff','Position',[870,130,120,25], 'Callback', @extractimg_click); %To extract image in roi to tiff
deltempbut = uicontrol('Style','pushbutton','String','Delete Temp Data','Position',[870,90,120,30], 'Callback', @deltemp_click); %Delete temporary data in database
closefigsbut = uicontrol('Style','pushbutton','String','Close sub figures','Position',[870,20,120,25], 'Callback', @closefig_click); %Close all figures except GUI
clearlastfitbut= uicontrol('Style','pushbutton','String','Clear last fit data','Position',[870,55,120,25], 'Callback', @clearlastfit_click); %Clears last fit data (to be used in case of change of cursor or angle)


%% Normalize the components
% Change units to normalized so components resize automatically.
f.Units = 'normalized';
img.Units = 'normalized';
updatebut.Units = 'normalized';
quickres.Units = 'normalized';
framelist.Units = 'normalized';
zon.Units = 'normalized';
zreset.Units = 'normalized';
curs.Units = 'normalized';
panbtn.Units = 'normalized';
rotatebtn.Units = 'normalized';
rotangleinput.Units = 'normalized';
colormapname.Units = 'normalized';
dblist.Units = 'normalized';
savebut.Units = 'normalized';
loadbut.Units = 'normalized';
delbut.Units = 'normalized';
add2anabut.Units = 'normalized';
nextimgname_text.Units = 'normalized';
nextimgname.Units = 'normalized';
imgid_text.Units = 'normalized';
imgidbox.Units = 'normalized';
imgname_text.Units = 'normalized';
imgname.Units = 'normalized';
timestmp_text.Units = 'normalized';
timestmpstatus.Units = 'normalized';
savestatus_text.Units = 'normalized';
savestatus.Units = 'normalized';
seqidstatus_text.Units = 'normalized';
seqidstatus.Units = 'normalized';
seqnamestatus_text.Units = 'normalized';
seqnamestatus.Units = 'normalized';
seqdescstatus_text.Units = 'normalized';
seqdescstatus.Units = 'normalized';
fitplt.Units = 'normalized';
fitbut.Units = 'normalized';
singlefitbut.Units = 'normalized';
cftoolbut.Units = 'normalized';
fitplotcheck.Units = 'normalized';
fitres.Units = 'normalized';
analysisimgidbox.Units = 'normalized';
showanalysisimgidbut.Units = 'normalized';
xvariable_text.Units = 'normalized';
xvarmode.Units = 'normalized';
xvarnumbermode.Units = 'normalized';
xvarnamemode.Units = 'normalized';
fitxvar.Units = 'normalized';
xvardbmode.Units = 'normalized';
xvardropmenu.Units = 'normalized';
yvariable_text.Units = 'normalized';
analysisdblist.Units = 'normalized';
load2anabut.Units = 'normalized';
del2anabut.Units = 'normalized';
clear2anabut.Units = 'normalized';
fitlist_text.Units = 'normalized';
fitlist.Units = 'normalized';
fittype_text.Units = 'normalized';
fittypelist.Units = 'normalized';
extractimgbut.Units = 'normalized';
deltempbut.Units = 'normalized';
closefigsbut.Units = 'normalized';
fitoutputnum_text.Units = 'normalized';
fitoutputnum.Units = 'normalized';
nextname_text.Units = 'normalized';
nextnameOption.Units = 'normalized';
clearlastfitbut.Units = 'normalized';
xvarincrement.Units = 'normalized';
xvardecrement.Units = 'normalized';
yvarincrement.Units = 'normalized';
yvardecrement.Units = 'normalized';
rotincrement.Units = 'normalized';
rotdecrement.Units = 'normalized';
autoupbox.Units = 'normalized';
addnext.Units = 'normalized';
shownamecheck.Units = 'normalized';
pcacbox.Units = 'normalized';
pcashowbox.Units = 'normalized';
imgModeText.Units = 'normalized';
imgMode.Units = 'normalized';
pcaTest.Units = 'normalized';


%% Initialize the UI
f.ToolBar = 'None'; %Hide the main toobar in GUI
f.MenuBar = 'None'; % Hide menu bar in GUI
f.CloseRequestFcn = @onClose;
zoom_img=zoom(img);
zoom_img.Enable = 'off';
% dbh=NET.Assembly('DatabaseHelper.dll');
updatefittypelist();
updatefitlist();
updatexvardropmenu();

%% Initialize Analysis
loadPCABasis();
pcaflag = false;


showimg(currentimgid);
while true 
    if ~ishandle(dblist)
        break;
    end
    updateimgidlist();
    updatedblist();
    pause(0.5);
    if (isvalid(pcashowbox) && get(pcashowbox,'Value') && ~pcaflag)
        showimg(currentimgid);
    end
end

if ishandle(dblist)         %To make sure no new figure opens on closing GUI, need to define global lines to move them around
    xl1 = line([xcurs(1) xcurs(1) ],[1 row]);
    xl2 = line([xcurs(2) xcurs(2)],[1 row]);
    yl1 = line([1 col], [ycurs(1) ycurs(1)]);
    yl2 = line([1 col], [ycurs(2) ycurs(2)]);
end



%% Defining all icons for toolbox for img
function [zoom_icon,pan_icon,curs_icon,rotate_icon] = icon_update()
    currentfolder = pwd;
    %Zoom icon
    [abc,~,alpha] = imread(fullfile(currentfolder,'tool_zoom_in.png'));
    zoom_icon = double(abc)/256/256;
    zoom_icon(~alpha) = NaN;
    %Pan icon
    [abc,~,alpha] = imread(fullfile(currentfolder,'tool_hand.png'));
    pan_icon = double(abc)/256/256;
    pan_icon(~alpha) = NaN;
    %Curs icon
    [abc,~,alpha] = imread(fullfile(currentfolder,'tool_data_cursor.png'));
    curs_icon = double(abc)/256/256;
    curs_icon(~alpha) = NaN;
    %Rotate icon
    [abc,~,alpha] = imread(fullfile(currentfolder,'tool_rotate_3d.png'));
    rotate_icon = double(abc)/256/256;
    rotate_icon(~alpha) = NaN;
end


%% Updating global imgidlist which is displayed in dblist
function updateimgidlist()
    m=max(cell2mat(imgidlist));
    sqlquery=['SELECT imageID FROM images ORDER BY imageID DESC LIMIT ' num2str(numImages)];
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    newimgidlist = curs1.Data;
    %If size is more than numImages then replace last with 1st in imgidlist
    length_imgidlist=length(imgidlist);
    if length_imgidlist > numImages
        for i=1:length(newimgidlist)
            imgidlist{i}=newimgidlist{i};
        end
    else
        imgidlist=newimgidlist;
    end  
    if max(cell2mat(imgidlist)) > m
        if get(nextnameOption,'Value') == 2
            nextname=get(nextimgname,'String');
            sqlquery2=['UPDATE images SET name ="', nextname,'" WHERE imageID = ', num2str(max(cell2mat(imgidlist)))];
            curs2=exec(conn, sqlquery2);
            close(curs2);
        elseif get(nextnameOption, 'Value') == 3
            namepart=get(nextimgname,'String');
            %Determine variable value to add to name string
            var_names=get(xvardropmenu,'String');
            var_value=get(xvardropmenu,'Value');
            var=var_names{var_value};            
            sqlquery3=['SELECT ',var,' FROM ciceroout WHERE runID = (SELECT runID_fk FROM images ORDER BY imageID DESC LIMIT 1)'];
            curs3=exec(conn, sqlquery3);
            curs3=fetch(curs3);
            varValue=num2str(cell2mat(curs3.Data));
            close(curs3);
            nextname=[namepart,' ',varValue];
            sqlquery2=['UPDATE images SET name ="', nextname,'" WHERE imageID = ', num2str(max(cell2mat(imgidlist)))];
            curs2=exec(conn, sqlquery2);
            close(curs2);
        end
        currentimgid=imgidlist(1);
        updatecurrentimginfo(currentimgid);
        showimg(currentimgid);
        [y,Fs] = audioread('sound.wav');
        sound(y,Fs);
        if get(addnext,'Value') == 1
            anaimgidlist=[anaimgidlist; currentimgid];
            updateanalysisdblist();
        end
        pcaflag = false;
    end
    close(curs1);
end    
    
%% Updates the fitting functions from specified folder 'funcnpath'
function updatefitlist()
    currentfolder = pwd;
    funcnpath=[currentfolder '\FittingFunctions'];
    funcns = dir([funcnpath '/*.m']);
    funcns = cellfun(@(x) x(1:end-2),{funcns.name},'Un',0);
    set(fitlist, 'string', funcns);
end


%% Loads 'pwoaSave.mat' to create the starting basis for PCA.
% Assumes that pwoSave is sorted by imageId ascending. I.E. That
% pwoaSave(1) is the oldest image.
function loadPCABasis()
    try 
        pwoa = 0;
        load('pwoaSave','pwoa');
        savedSize = size(pwoa,1);
        tempList = reshape(pwoa,savedSize,col*row);
        pwoa = tempList(max(savedSize-pcaBasisSize+1,1):savedSize,:);
    catch
        warning('No PCA basis found. Starting with 0.');
    end
    try
        load('oldestImgSave','oldestImgIndex');
    catch
        warning('No oldestImgSave.m found. Starting with 1');
    end
    
    eigenBasis = makeEigenBasis(pwoa);
end

%% Makes an eigenbasis of the images in the argument.
function [eigenOut] = makeEigenBasis(imageList)
    if size(imageList,1) < 2
        eigenOut = zeros(col*row,1);
    else
        eigenOut = pca(pwoa);
    end
end

%% Update analysis database list
function updateanalysisdblist()
    anaimgidlist=num2cell(unique(cell2mat(anaimgidlist),'stable'));
    anaimgidlist=num2cell(sort(cell2mat(anaimgidlist), 'descend'));
    sqlquery = ['SELECT name FROM images WHERE imageID IN (', strjoin(cellstr(num2str(cell2mat(anaimgidlist))),','),') ORDER BY imageID DESC'];
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    db = curs1.Data;
    close(curs1);
    set(analysisdblist, 'string', db);
    set(analysisdblist, 'Value', []);
end

%% Updates database list
function updatedblist()
    % Sql code for getting name of images
    if get(shownamecheck,'Value') == 1
        sqlquery = ['SELECT sequence.name, images.name FROM sequence, images WHERE (sequence.sequenceID = images.sequenceID_fk AND images.imageID IN(', strjoin(cellstr(num2str(cell2mat(imgidlist))),','),')) ORDER BY images.imageID DESC'];
        curs1=exec(conn, sqlquery);
        curs1=fetch(curs1);
        db = curs1.Data;
        close(curs1);
        db = strcat(db(:,1),{': '},db(:,2));
    else
        sqlquery = ['SELECT name FROM images WHERE imageID IN (', strjoin(cellstr(num2str(cell2mat(imgidlist))),','),') ORDER BY imageID DESC'];
        curs1=exec(conn, sqlquery);
        curs1=fetch(curs1);
        db = curs1.Data;
        close(curs1);
    end
    
    if ~isempty(anaimgidlist)
        for i = 1:length(db)
            if any(cell2mat(imgidlist(i))==cell2mat(anaimgidlist))
                db(i) = strcat('*',db(i));
            end
        end
    end
    set(dblist, 'string', db);
end

%% Call back function for clicking in database list
function dblist_click(~, ~)
    sel=get(dblist, 'Value');
    l_current=length(sel);
    l_past=length(selecimgidlist);
    if isempty(selecimgidlist) == 1
        selecimgidlist=imgidlist(sel);
        currentimgid=imgidlist(sel);
    elseif l_current == 1
        selecimgidlist = imgidlist(sel);
        currentimgid=imgidlist(sel);
    elseif l_current-l_past == 1
        currentimgid=num2cell(setxor(cell2mat(selecimgidlist), cell2mat(imgidlist(sel))));
        if length(currentimgid) > 1
            currentimgid=imgidlist(sel);
        end
        selecimgidlist=imgidlist(sel);
    else
        selecimgidlist = imgidlist(sel);
        currentimgid=imgidlist(sel);
    end
    if length(currentimgid) == 1
        updatecurrentimginfo(currentimgid);
        showimg(currentimgid);
    end
end    

%% Call back function for pressing key in database list
function dblist_keypress(source, eventdata)
    key = get(gcf,'CurrentKey');
    if(strcmp (key , 'delete'))
        del_click(source, eventdata);
    elseif(strcmp (key , 'rightarrow'))
        add2ana_click(source, eventdata);
    end
end    

%% Function to update text boxes showing details of currently/last selected image from dblist
function updatecurrentimginfo(currentimgid)
%ID
    set(imgidbox, 'string', num2str(cell2mat(currentimgid)));
%Name
    sqlquery=['SELECT name FROM images WHERE imageID =', num2str(cell2mat(currentimgid))];
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    info = curs1.Data;
    close(curs1);
    set(imgname, 'string', info);
%Timestamp
    sqlquery=['SELECT timestamp FROM images WHERE imageID =', num2str(cell2mat(currentimgid))];
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    info1 = curs1.Data;
    close(curs1);
    set(timestmpstatus, 'string', info1);
%Save
    sqlquery=['SELECT type FROM images WHERE imageID =', num2str(cell2mat(currentimgid))];
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    info = curs1.Data;
    close(curs1);
    if strcmp(info,'perm')
        set(savestatus, 'string', 'Yes');
    else
        set(savestatus, 'string', 'No');
    end
%SeqID
    sqlquery=['SELECT sequenceID_fk FROM images WHERE imageID =', num2str(cell2mat(currentimgid))];
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    seqID = curs1.Data;
    close(curs1);
    set(seqidstatus, 'string', seqID);
%SeqName
    sqlquery=['SELECT name FROM sequence WHERE sequenceID =', num2str(cell2mat(seqID))];
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    info2 = curs1.Data;
    close(curs1);
    set(seqnamestatus, 'string', info2);
%SeqDesc
    sqlquery=['SELECT description FROM sequence WHERE sequenceID =', num2str(cell2mat(seqID))];
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    info3 = curs1.Data;
    close(curs1);
    set(seqdescstatus, 'string', info3);
end


%% Call back function for clicking in frame list
function framelist_click(~, ~)
    showimg(currentimgid);    
end

%% Call back function for clicking PCA checkbox
function pca_click(~, ~)
    showimg(currentimgid);    
end

%% Call back function for change in capture mode 'hdwmode'
function hdwmode_change (~, ~)
    showimg(currentimgid);
end

%% Call back function for change in colormap mode
function colormapname_click (~, ~)
    showimg(currentimgid);
end

%% Gets the image mode requested. 1 = normal, 2 = kinetics with two frames, 3 = kinetics all on one frame.
function imgmodeout = getImgMode()
    imgmodeout = imgMode.Value;
end

%% For updating image in 'img', inserted one input parameter as it has to be globally defined
function showimg(filenum)
    imgid=cell2mat(filenum(1));
    framenum=get(framelist, 'Value');
    imgmode = getImgMode();
    
    r = getImage(imgid,imgmode,framenum);
    cla(img);
    [~]= imagesc(r, 'Parent', img);
    if get(colormapname,'Value') == 1
        load('MyColormaps','mycmap')
        colormap(img, mycmap);
    else
        colormap(img, gray);
    end
    if framenum == 1
        caxis(img, [0 1.2]);
        [~] = colorbar(img,'XTickLabel',{'0' '1.20'}, 'XTick', [0 1.2]);
    else
        caxis(img, [min(r(:)) max(r(:))]);
        [~] = colorbar(img,'XTickLabel',{num2str(min(r(:))) num2str(max(r(:)))},'XTick', [min(r(:)) max(r(:))]);
    end

    curs_update();
    
    if get(autoupbox,'Value') == 1
        update_but(0,0);
    end
end

%% Performs PCA on the selected PWA image using the global eigenBasis.
function [imgOut] = doPCA(imgIn)
    
    toImg = reshape(imgIn,size(imgIn,1)*size(imgIn,2),1);
    if length(toImg) < col*row
        toImg(col*row,1) = 0; % add padding
    end
    
    meanImg = mean(pwoa,1)';

    c = (toImg-meanImg)'*eigenBasis;
    estPWOA = (eigenBasis*c'+meanImg);
        
    toImg = toImg(1:size(imgIn,1)*size(imgIn,2));
    estPWOA = estPWOA(1:size(imgIn,1)*size(imgIn,2));
    
    imgOut = reshape(toImg./estPWOA,size(imgIn,1),size(imgIn,2));

end

%% Adds a new col*row image to the PCA basis and updates the eigenbasis.
function addToBasis(newImg)
    X = reshape(newImg,1,size(newImg,1)*size(newImg,2));
    if length(X) < col*row
        X(1,col*row) = 0;
    end
    pwoa(oldestImgIndex,:) = X;
    oldestImgIndex = mod(oldestImgIndex, pcaBasisSize) + 1;
    eigenBasis = makeEigenBasis(pwoa);
end

%% Gets the requested image from the database.
function [r] = getImage(imgid,imgmode,framenum)
    if ((get(pcacbox,'Value') || get(pcashowbox,'Value')) && framenum == 1) % PCA is checked and we want the absorption image
        sqlquery2=['SELECT pcadata FROM images WHERE imageID = ', num2str(imgid)];
        curs2=exec(conn, sqlquery2);
        curs2=fetch(curs2);
        bdata=curs2.Data;
        close(curs2); 
    end
    
    if pcaTest.Value
        bdata = mat2cell('null',1);
    end
    
    if (get(pcacbox,'Value') && framenum == 1 && strcmp('null',cell2mat(bdata))) % We haven't done PCA on this image yet.
        sqlquery2=['SELECT data, cameraID_fk FROM images WHERE imageID = ', num2str(imgid)];
        curs2=exec(conn, sqlquery2);
        curs2=fetch(curs2);
        bdata=curs2.Data;
        close(curs2);
        blobdata=typecast(cell2mat(bdata(1)),'int16');
        sqlquery3=['SELECT cameraHeight, cameraWidth, Depth FROM cameras WHERE cameraID = ', num2str(cell2mat(bdata(2)))];
        curs3=exec(conn, sqlquery3);
        curs3=fetch(curs3);
        camdata=curs3.Data;
        close(curs3);
        camdata=cell2mat(camdata);
        s=[camdata(1),camdata(2),camdata(3)];            
        a=Blob2Matlab(blobdata,s);

        newPWOA = data_evaluation(a, imgmode, 3);
        newPWA = data_evaluation(a, imgmode, 2);
        newDark = data_evaluation(a, imgmode, 4);

        newPWA = double(max(min(newPWA - newDark,65535),1));
        newPWOA = double(max(min(newPWOA - newDark,65535),1));

        addToBasis(newPWOA);
        r = doPCA(newPWA);

        if ~pcaTest.Value
        tableName = 'images';
        colName = {'pcadata'};
        data = {typecast(reshape(r,1,col*row),'int8')};
        whereClause = ['WHERE imageID = ', num2str(imgid)];
        update(conn,tableName,colName,data,whereClause);
        end
        pcaflag = true;
        
    elseif ((get(pcacbox,'Value') || get(pcashowbox,'Value')) && framenum == 1 && ~strcmp('null',cell2mat(bdata)))
        blobdata=typecast(cell2mat(bdata),'double');
        s=[col row];
        r=Blob2Matlab(blobdata,s);
        
        pcaflag = true;
    
    else
        % normal imaging
        sqlquery2=['SELECT data, cameraID_fk FROM images WHERE imageID = ', num2str(imgid)];
        curs2=exec(conn, sqlquery2);
        curs2=fetch(curs2);
        bdata=curs2.Data;
        close(curs2);
        blobdata=typecast(cell2mat(bdata(1)),'int16');
        sqlquery3=['SELECT cameraHeight, cameraWidth, Depth FROM cameras WHERE cameraID = ', num2str(cell2mat(bdata(2)))];
        curs3=exec(conn, sqlquery3);
        curs3=fetch(curs3);
        camdata=curs3.Data;
        close(curs3);
        camdata=cell2mat(camdata);
        s=[camdata(1),camdata(2),camdata(3)];            
        a=Blob2Matlab(blobdata,s);

        r = data_evaluation(a, imgmode,framenum);
        
    end

    % rotate
    button_state = get(rotatebtn,'Value');
    if button_state == get(rotatebtn,'Max')
        r=imrotate(r, str2double(rotangle), 'bilinear', 'loose');          %Imrotate type and crop property can be defined here
    end
end


%% Gets the ROI as set by the cursors.
function [roi] = get_roi(imgid,imgmode,framenum)
    framedata=getImage(imgid,imgmode, framenum);
    minx=round(min(xcurs));
    maxx=round(max(xcurs));
    miny=round(min(ycurs));
    maxy=round(max(ycurs));
    roi=framedata(miny:maxy,minx:maxx);           %Swapped x & y due to matlab notation

end

%% Call back function for load button for dblist. Can only load 1 imageID now
function load_click(~, ~)
    newid=inputdlg('Enter Image ID:', 'Load');
    a=cell2mat(newid);
    l=length(imgidlist);
    imgidlist(l+1)=num2cell(str2double(a));
    imgidlist=num2cell(unique(cell2mat(imgidlist),'stable'));
    updatedblist();
end

%% Call back function for save button for dblist. Change type from temp to perm
function save_click(~, ~)
    sel=get(dblist, 'Value');
    id= imgidlist(sel);
    sqlquery2=['UPDATE images SET type = ''perm'' WHERE imageID IN (', strjoin(cellstr(num2str(cell2mat(id))),','),')'];
    curs2=exec(conn, sqlquery2);
    close(curs2);
    updatecurrentimginfo(currentimgid(1));
end    

%% Call back function for delete button for dblist
function del_click(~, ~)
    sel=get(dblist, 'Value');
    id= imgidlist(sel);
    sqlquery2 = ['DELETE FROM images WHERE imageID IN (', strjoin(cellstr(num2str(cell2mat(id))),','),')'];
    curs2=exec(conn, sqlquery2);
    close(curs2);
    imgidlist(sel)=[];
    updateimgidlist();
    updatedblist();
    set(dblist, 'Value', 1);
    updatecurrentimginfo(currentimgid);
    for i=1:length(anaimgidlist)        %Removes deleted id from anaimgidlist if it is there
        if anaimgidlist{i} == cell2mat(id)
            anaimgidlist(i)=[];
        end
    end
    updateanalysisdblist();
end    

%% Function for renaming image in dblist
function rename_imgname(~, ~)
    sel=get(dblist,'Value');
    newname=get(imgname,'String');
    id= cell2mat(imgidlist(sel));
    sqlquery2=['UPDATE images SET name ="', num2str(cell2mat(newname)),'" WHERE imageID = ', num2str(id)];
    curs2=exec(conn, sqlquery2);
    close(curs2);
    updateimgidlist();
    updatedblist();
    updatecurrentimginfo(currentimgid);
end

%% Call back function for add 2 analysis button for dblist. 
function add2ana_click(~, ~)
    sel=get(dblist,'Value');
    anaimgidlist=[anaimgidlist; imgidlist(sel)];
    updateanalysisdblist();
end

%% Call back function for clicking in database list
function analysisdblist_click(~, ~)
    sel=get(analysisdblist, 'Value');
    id=cell2mat(anaimgidlist(sel));
    set(analysisimgidbox,'String',num2str(id));
end

%% Call back function for pressing key in analysisdb list
function analist_keypress(source, eventdata)
    key = get(gcf,'CurrentKey');
    if(strcmp (key , 'delete'))
        del2ana_click(source, eventdata);
    end
end

%% Call back function for load button for analysisdblist
function load2ana_click(~, ~)
    newid=inputdlg('Enter Sequence ID:', 'Load images for analysis');
    a=cell2mat(newid);
    sqlquery2=['SELECT imageID FROM images WHERE sequenceID_fk = ', a];
    curs2=exec(conn, sqlquery2);
    curs2=fetch(curs2);
    imgids=curs2.Data;
    close(curs2);
    anaimgidlist = [anaimgidlist; imgids];
    updateanalysisdblist();
end

%% Call back function for delete from analysisdblist. 
function del2ana_click(~, ~)
    selected = get(analysisdblist,'Value');
%     list = get(analysisdblist, 'String');
    anaimgidlist(selected)=[];
    updateanalysisdblist();
end

%% Call back function for clearing all images from analysisdblist. 
function clear2ana_click(~, ~)
    anaimgidlist = [];
    updateanalysisdblist();
end

%% Zoom button (currently handles all zoom requirements)
function zoom_on (~,~)
    button_state = get(zon,'Value');
    if button_state == get(zon,'Max')
        zoom_img.Enable = 'on';
    elseif button_state == get(zon,'Min')
        zoom_img.Enable = 'off';
    end
end

%Zoom reset to 100%
function zoom_reset(~, ~)
    axis(img, [1 col 1 row]);
    set(zreset, 'Value', 0);
end

%% Turn pan on for the img axes
function pan_on(~, ~)
    button_state = get(panbtn,'Value');
    if button_state == get(panbtn,'Max')
        pan(img,'On');
    elseif button_state == get(panbtn,'Min')
        pan(img,'Off');
    end
end

%% Rotate for img axes, we just do showimg() with rotation
function rotate_on(~, ~)
    showimg(currentimgid);
end

%% Callback function to enter new angle of rotation for img axes
function enter_rotangle(~, ~)
    rotangle=get(rotangleinput,'String');
    showimg(currentimgid);
end

%% Callback function to increment rotangle
function rotinc_click(~,~)
    temp = get(rotangleinput,'String');
    rotangle = num2str(str2double(temp) + 1);
    set(rotangleinput,'String',rotangle);
    showimg(currentimgid);
end

%% Callback function to decrement rotangle
function rotdec_click(~,~)
    temp = get(rotangleinput,'String');
    rotangle = num2str(str2double(temp) - 1);
    set(rotangleinput,'String',rotangle);
    showimg(currentimgid);
end


%% Cursor button call back function 
function curs_on (~, ~) 
    button_state = get(curs,'Value');
    if button_state == get(curs,'Max')
        [xcurs,ycurs] = ginputax(img,2);    
        showimg(currentimgid);
        set(curs,'Value',0);
    end
end

%% Cursor update function
function curs_update()
    axes(img)
    hold(img,'on');
        xl1 = line([xcurs(1) xcurs(1) ],[1 row]);
        xl2 = line([xcurs(2) xcurs(2)],[1 row]);
        yl1 = line([1 col], [ycurs(1) ycurs(1)]);
        yl2 = line([1 col], [ycurs(2) ycurs(2)]);
        draggable(xl1,'constraint','horizontal', 'endfcn', @xl1_change);
        draggable(xl2,'constraint','horizontal', 'endfcn', @xl2_change);
        draggable(yl1,'constraint','vertical', 'endfcn', @yl1_change);
        draggable(yl2,'constraint','vertical', 'endfcn', @yl2_change);
end

function xl1_change(~, ~)
    x1=get(xl1, 'XData');
    x3=get(xl2, 'XData');
    if x3(1) == xcurs(1)
        xcurs(2) = x1(1);
    else
        xcurs(1) = x1(1);
    end
end
    
function xl2_change(~, ~)
    x1=get(xl1, 'XData');
    x3=get(xl2, 'XData');
    if x1(1) == xcurs(1)
        xcurs(2) = x3(1);
    else
        xcurs(1) = x3(1);
    end
end

function yl1_change(~, ~)
    y1=get(yl1, 'YData');
    y3=get(yl2, 'YData');
    if y3(1) == ycurs(1)
        ycurs(2) = y1(1);
    else
        ycurs(1) = y1(1);
    end
end

function yl2_change(~, ~)
    y1=get(yl1, 'YData');
    y3=get(yl2, 'YData');
    if y1(1) == ycurs(1)
        ycurs(2) = y3(1);
    else
        ycurs(1) = y3(1);
    end
end


%% Callback function for update button for quick results
function update_but(~, ~)
    framenum=get(framelist, 'Value');
    imgmode = getImgMode();
    imgid=cell2mat(currentimgid);
    data_roi = get_roi(imgid,imgmode,framenum);
   
    try
        currentfolder = pwd;
        fitpath=[currentfolder '\FittingFunctions'];
        addpath(fitpath);
        norm_n = NormN_Count(data_roi,0);
        xGauss = Gauss_Fit(data_roi,0);
        yGauss = Gauss_Fit(data_roi',0);
    catch
        warning('Error with fitting.');
        norm_n = 0;
        xGauss = [0,0,0];
        yGauss = [0,0,0];
    end
    set(quickres, 'String', {['Norm N Count: ' num2str(norm_n)]; ['X Width: ' num2str(2^(3/2) * xGauss(3))]; ['Y Width: ' num2str(2^(3/2) * yGauss(3))]; ['X Center: ' num2str(xGauss(2))]; ['Y Center: ' num2str(yGauss(2))]; ['X Curs: ' num2str(round(min(xcurs))) ' , ' num2str(round(max(xcurs)))];['Y Curs: ' num2str(round(min(ycurs))) ' , ' num2str(round(max(ycurs)))];});
end

%% Callback function for fit button to plot results in fitplt
function fit_click(~, ~)
    fitnum=get(fitlist,'Value');
    currentfolder = pwd;
    fitpath=[currentfolder '\FittingFunctions'];
    addpath(fitpath);
    fbase=dir([fitpath '/*.m']);
    fitfuncname=fbase(fitnum);
    fitfunc=str2func(fitfuncname.name(1:end-2));

    framenum=get(framelist, 'Value');
    imgmode = getImgMode();
    if get(fitplotcheck, 'Value') == get(fitplotcheck, 'Max')   %Checks whether to plot sub-figure for each fit
        eachplot = 1;
    else
        eachplot = 0;
    end   
    random_output=fitfunc(rand(50),0);
    resy=zeros(length(anaimgidlist),length(random_output));
    %To reduce redundancy by avoiding calculation of fitting function for images fitted in last run
    temp_analysislist=anaimgidlist;
    last_ids=getappdata(f,'ids');
    last_resy=getappdata(f,'resulty');
    for i=1:length(temp_analysislist)
        for j=1:length(last_ids)
            if temp_analysislist{i}== last_ids{j}
                resy(i,:)=last_resy(j,:);
                temp_analysislist{i}=0;
            end
        end
    end    
    %Calculation of y variable
    for i=1:length(temp_analysislist)
        if temp_analysislist{i} ~= 0 
            imgid=cell2mat(anaimgidlist(i));
            data_roi = get_roi(imgid,imgmode,framenum);
            resy(i,:)=fitfunc(data_roi, eachplot);             %Input is absorption image, not Optical Density
        end
    end    
    %Calculation of x variable
    resx=xvardet(length(anaimgidlist));    %Determines value/array of x variable for fit plot
    % Fittype inclusion
    fittyp=getfittype();
    output_num=int32(str2double(get(fitoutputnum, 'String')));
    if strcmp(fittyp,'') == 1
        axes(fitplt);
        resulty=resy(:,output_num)';
        scatter(resx, resulty);
        showdata1='';
        for i=1:length(resulty)
            showdata1{i}=sprintf('%s , %s', num2str(resx(i)), num2str(resulty(i)));
        end
        showdata=sprintf('\n%s', showdata1{:});
        set(fitres,'String',{'Results:';['Data: ' showdata]});
    elseif strcmp(fittyp, 'bessel2') == 1       %Bessel2 just works for Kapitza Dirac as of now
        leastsquare=cell(length(anaimgidlist),1);
        ls=@(theta) 0;
        for i=1:length(anaimgidlist)
            leastsquare{i}=@(theta)(besselj(0,theta*resx(i))^2-resy(i,1))^2+(besselj(1,theta*resx(i))^2-resy(i,2))^2+(besselj(2,theta*resx(i))^2-resy(i,3))^2;
            ls=@(theta) ls(theta)+leastsquare{i}(theta);
        end
        theta0=0.2;
        theta=fminsearch(ls, theta0);
        volts_5recoil=(0.0000125*1000*2*pi*5*2.02781)/(2*theta);
        resultx = zeros(1,length(anaimgidlist));
        resulty0 = zeros(1,length(anaimgidlist));
        resulty1 = zeros(1,length(anaimgidlist));
        resulty2 = zeros(1,length(anaimgidlist));
        for i=1:length(anaimgidlist)
            resulty0(i)=resy(i,1);
            resulty1(i)=resy(i,2);
            resulty2(i)=resy(i,3);
            resultx(i)=theta*resx(i);
        end
        showdata1='';
        for i=1:length(resulty0)
            showdata1{i}=sprintf('%s , %s', num2str(resx(i)), num2str(resulty0(i)));
        end
        showdata=sprintf('\n%s', showdata1{:});
        axes(fitplt);
        plot(resultx,resulty0,'bo',resultx,resulty1,'ro',resultx,resulty2,'go');
        hold on
        fplot(@(x) besselj(0,x)^2, [min(resultx) max(resultx)],'b');
        fplot(@(x) besselj(1,x)^2, [min(resultx) max(resultx)],'r');
        fplot(@(x) besselj(2,x)^2, [min(resultx) max(resultx)],'g');
        hold off
        set(fitres,'String',{'Results:';['Theta: ' num2str(theta)];['Volts/5 Erec: ' num2str(volts_5recoil)];['Bessel 0 Data: ' showdata]});
    else
        resulty=resy(:,output_num)';
        fit_data=fit(resx',resulty',fittyp);
        axes(fitplt);
        plot(fit_data, resx, resulty,'o');
        coeff_name=coeffnames(fit_data);
        coeff_value=coeffvalues(fit_data);
        coeff1='';
        for i=1:length(coeff_name)
            coeff1{i}=sprintf('%s = %s', coeff_name{i}, num2str(coeff_value(i)));
        end
        coeff=sprintf('\n%s', coeff1{:});
        showdata1='';
        for i=1:length(resulty)
            showdata1{i}=sprintf('%s , %s', num2str(resx(i)), num2str(resulty(i)));
        end
        showdata=sprintf('\n%s', showdata1{:});
        set(fitres,'String',{'Results:';['Formula: ' formula(fit_data)]; ['Coefficients: ' coeff]; ['Data: ' showdata]});
    end
    setappdata(f, 'ids', anaimgidlist);
    setappdata(f, 'resulty', resy);
    setappdata(f, 'resultx', resx);
end

%% Callback function for fit button to plot results in fitplt
function singlefit_click(~, ~)
    fitnum=get(fitlist,'Value');
    currentfolder = pwd;
    fitpath=[currentfolder '\FittingFunctions'];
    addpath(fitpath);
    fbase=dir([fitpath '/*.m']);
    fitfuncname=fbase(fitnum);
    fitfunc=str2func(fitfuncname.name(1:end-2));

    framenum=get(framelist, 'Value');
    imgmode = getImgMode();
    eachplot = 1;
    val= get(analysisdblist,'Value');
    [~]= get(analysisdblist,'String');
    imgid=cell2mat(anaimgidlist(val));
    p = get_roi(imgid,imgmode,framenum);
    resy=fitfunc(p, eachplot);
    
    output_num=int32(str2double(get(fitoutputnum, 'String')));
    resulty=resy(:,output_num)';
    % X Variable for selected imageID 
    xvarmodesel=get(xvarmode, 'SelectedObject');
    xvarmodevalue=xvarmodesel.Tag;
    resx=zeros(1);
    if xvarmodevalue == '1'
        resx=1;
    elseif xvarmodevalue == '2'
        filenames=get(analysisdblist, 'String');
        i_sel=get(analysisdblist, 'Value');
        name=filenames{i_sel};
        splitfilename=strsplit(name);
        b=cellfun(@str2double,splitfilename(:),'un',0).';
        v=int32(str2double(xvar));
        na=cell2mat(b(min(v,length(b))));
        if isnan(na)
            na = 0;
        end
        resx=na;
    elseif xvarmodevalue == '3'
        i_sel=get(analysisdblist, 'Value');
        imgid=cell2mat(anaimgidlist(i_sel));
        var_names=get(xvardropmenu,'String');
        var_value=get(xvardropmenu,'Value');
        var=var_names{var_value};            
        sqlquery3=['SELECT ',var,' FROM ciceroout WHERE runID = (SELECT runID_fk FROM images WHERE imageID = ', num2str(imgid), ')'];
        curs3=exec(conn, sqlquery3);
        curs3=fetch(curs3);
        resx=cell2mat(curs3.Data);
        close(curs3);
    end
    
    showdata=sprintf('\n%s %s', num2str(resx), num2str(resulty));
    set(fitres,'String',{'Results:';['Data: ' showdata]});
end

%% Determines value/array of x variable for fit plot
function resx = xvardet(n)
    xvarmodesel=get(xvarmode, 'SelectedObject');
    xvarmodevalue=xvarmodesel.Tag;
    resx=zeros(1,n);
    if xvarmodevalue == '1'
        for i=1:n
            resx(i)=i;
        end
    elseif xvarmodevalue == '2'
        %From file name (can also get name from querying SQL but this should be faster)
        for i=1:n
            filenames=get(analysisdblist, 'String');
            name=filenames{i};
            splitfilename=strsplit(name);
            b=cellfun(@str2double,splitfilename(:),'un',0).';
            v=int32(str2double(xvar));
            na=cell2mat(b(min(v,length(b))));
            if isnan(na)
                na = 0;
            end
            resx(i)=na;
        end
    elseif xvarmodevalue == '3'
        %From DB (currently it only gets variable from ciceroout, need case statement to include all tables)
        for i=1:n
            imgid=cell2mat(anaimgidlist(i));
            var_names=get(xvardropmenu,'String');
            var_value=get(xvardropmenu,'Value');
            var=var_names{var_value};            
            sqlquery3=['SELECT ',var,' FROM ciceroout WHERE runID = (SELECT runID_fk FROM images WHERE imageID = ', num2str(imgid), ')'];
            curs3=exec(conn, sqlquery3);
            curs3=fetch(curs3);
            resx(i)=cell2mat(curs3.Data);
            close(curs3);
        end
    end
end

%% Opens cftool with data from results
function cftool_click(~, ~)
    output_num=int32(str2double(get(fitoutputnum, 'String')));
    total_y=getappdata(f,'resulty');
    x=getappdata(f,'resultx');
    y=total_y(:, output_num);
    cftool(x,y);
end


%% Callback function to enter position of x variable in name
function enter_fitxvar(~, ~)
    xvar=get(fitxvar,'String');
end

%% Callback function to increment xvar
function xvarinc_click(~,~)
    tempx = get(fitxvar,'String');
    xvar = num2str(str2double(tempx) + 1);
    set(fitxvar,'String',xvar);
end

%% Callback function to decrement xvar
function xvardec_click(~,~)
    tempx = get(fitxvar,'String');
    xvar = num2str(max(str2double(tempx) - 1,1));
    set(fitxvar,'String',xvar);
end


%% Callback function to increment yvar
function yvarinc_click(~,~)
    tempy = get(fitoutputnum,'String');
    yvar = num2str(str2double(tempy) + 1);
    set(fitoutputnum,'String',yvar);
end

%% Callback function to decrement yvar
function yvardec_click(~,~)
    tempy = get(fitoutputnum,'String');
    yvar = num2str(max(str2double(tempy) - 1,1));
    set(fitoutputnum,'String',yvar);
end


%% To update list of all variables in database (SQL can be edited to make it a particular table instead)
function updatexvardropmenu()
    sqlquery='SELECT COLUMN_NAME FROM information_schema.columns WHERE table_name="ciceroOut" AND COLUMN_NAME NOT LIKE "_unusedcolumn_%"';
    curs1=exec(conn, sqlquery);
    curs1=fetch(curs1);
    columndata = curs1.Data;
    close(curs1);
    set(xvardropmenu, 'String', columndata);
end

%% 
function updatefittypelist()
    set(fittypelist, 'String', {'Gaussian','Exponential','Sine','Line','Bessel2','None'});
end

%% 
function [type] = getfittype()
    a=get(fittypelist,'Value');
    switch a
            case 1
                type='gauss1';
            case 2
                type='exp1';
            case 3
                type='sin1';
            case 4
                type='poly1';
            case 5
                type='bessel2';
            case 6
                type='';  %No fit type case. Add custom fit type here & updatefittypelist()
     end
end

%% Extract specified ROI of image shown in analyzer to tiff, opens prompt to specify destination
function extractimg_click(~, ~)
    framenum=get(framelist, 'Value');
    imgmode = getImgMode();
    imgid=cell2mat(currentimgid);
    data = get_roi(imgid,imgmode,framenum);
    
    [filename, pathname] = uiputfile('.tiff','Save as');

    imwrite(data,'temp.png','png'); %Need to write it as png and transfer from png to tiff since writing directly to tiff throws error 
    c=imread('temp.png');
    imwrite(c,[pathname '\' filename],'tiff');
end

%% Delete all imageID marked as Temp (should be done every few days)
function deltemp_click(~, ~)
    choice=questdlg('Are you sure you want to delete temporary images?', 'ATTENTION', 'Yes!', 'No', 'No');
    switch choice
        case 'Yes'
            sqlquery2='DELETE FROM images WHERE type="temp"';
            curs2=exec(conn, sqlquery2);
            close(curs2);
    end            
end

%% Close all figures except GUI
function closefig_click(~, ~)
    figs = get(0,'children');
    figs(figs == gcf) = []; % delete your current figure from the list
    close(figs)    
end

%% Show selected image in analysisdblist on img axes
function showanalysisimgid_click(~, ~)
    val= get(analysisdblist,'Value');
    currentimgid=anaimgidlist(val);
    showimg(currentimgid);
end

%% Clears data saved from last fit (required when fit after changing cursor or angle)
function clearlastfit_click(~, ~)
    setappdata(f, 'ids', []);
end

%% To update string of nextimgname edit box on pressing enter 
function nextimgname_enter(~, eventdata)
    k=eventdata.Key;
    if strcmp(k,'return')
        uicontrol(nextimgname_text);        %Changes focus to other object than edit textbox to update nexxtimgname string
    end
end

%% To update string of imgname edit box on pressing enter 
function imgname_enter(~, eventdata)
    k=eventdata.Key;
    if strcmp(k,'return')
        uicontrol(imgname_text);        %Changes focus to other object than edit textbox to update imgname string
    end
end

%% To maneuver with keyboard shortcuts
function globalShortcuts(source, eventdata)
    k=eventdata.Key;
    if strcmp(k,'f5')
        update_but(source,eventdata);
    elseif strcmp(k,'f9')
        fit_click(source,eventdata);
    elseif strcmp(k,'insert')
        uicontrol(nextimgname);
    end
end

%% On close, save settings.
function onClose(~, ~)
    if get(pcacbox,'Value') % PCA is checked
        save('pwoaSave','pwoa');
        save('oldestImgSave','oldestImgIndex');
    end
    delete(f);
end

end