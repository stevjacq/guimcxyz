 function GUI_mcxyz
% function gui_mcxyz.m
% Uses:
%   maketissue_test.m
%       makeTissueList_hf.m
%       makecmap_hf.m
%   edit_table.m
%
% updated 
%   6sep: fix GOMCXYZ
%   19sep: add skin temp
clear all; clc; home
gui = 'GUI_mcxyz';
disp(gui)

homedir = pwd;
cd(homedir)
dirname = homedir;
addpath(homedir)

%% default parameters
fid = [];
dlmtr = '/';


simname = 'pick_a_name';
tissuename = 'skin_vessel'; 

u = [];
rnames = '';
cnames = '';
SHOWflag = 0;
GOflag = 0;
Tbx = [];
Tez = [];
x = [];
y = [];
z = [];
Fyxz = [];
Qyxz = [];
Qtt = [];
Up = [];
sfig1 = '';
sfig2 = '';
clr(1,:) = [1 1 1]*0.99; % white
clr(2,:) = [.7 1 .7]; % light green
clr(3,:) = [1 1 1]*0.7; % gray
clr(4,:) = [1 0 0]; % red
nm = [];

i=1;
mc(i).name = sprintf('%d  time_min',i); 
mc(i).param = 0.25;
i=i+1;
mc(i).name = sprintf('%d  Nx',i); 
mc(i).param = 200;
i=i+1;
mc(i).name = sprintf('%d  Ny',i); 
mc(i).param = 200;
i=i+1;
mc(i).name = sprintf('%d  Nz',i); 
mc(i).param = 200;
i=i+1;
mc(i).name = sprintf('%d  dx',i); 
mc(i).param = 0.0100;
i=i+1;
mc(i).name = sprintf('%d  dy',i); 
mc(i).param = 0.0100;
i=i+1;
mc(i).name = sprintf('%d  dz',i); 
mc(i).param = 0.0100;
i=i+1;
mc(i).name = sprintf('%d  mcflag',i); 
mc(i).param = 0;
i=i+1;
mc(i).name = sprintf('%d  boundaryflag',i); 
mc(i).param = 2;
i=i+1;
mc(i).name = sprintf('%d  xs',i); 
mc(i).param = 0;
i=i+1;
mc(i).name = sprintf('%d  ys',i); 
mc(i).param = 0;
i=i+1;
mc(i).name = sprintf('%d  zs',i); 
mc(i).param = 0.250;
i=i+1;
mc(i).name = sprintf('%d  xfocus',i); 
mc(i).param = 0;
i=i+1;
mc(i).name = sprintf('%d  yfocus',i); 
mc(i).param = 0;
i=i+1;
mc(i).name = sprintf('%d  zfocus',i); 
mc(i).param = 0.0500;
i=i+1;
mc(i).name = sprintf('%d  radius',i); 
mc(i).param = 0.5;
i=i+1;
mc(i).name = sprintf('%d  waist',i); 
mc(i).param = 0.010;
i=i+1;
mc(i).name = sprintf('%d  zsurf',i); 
mc(i).param = 0.25;
i=i+1;
mc(i).name = sprintf('%d  launchflag',i); 
mc(i).param = 0;
i=i+1;
mc(i).name = sprintf('%d  ux0',i);
mc(i).param = 0.25;
i=i+1;
mc(i).name = sprintf('%d  uy0',i);
mc(i).param = 0.25;
i=i+1;
mc(i).name = sprintf('%d  yz0',i); 
mc(i).param = 0.25;

%%%%%%%%%%%
% control panel
%% %
sz0 = get(0,'screensize');
sw0 = sz0(3);
sh0 = sz0(4);
whr = sw0/sh0; % screen width / height ratio

%%% set figure windows
% Figure 1
scale = [1 .83]*0.5;     % horiz/vert size
pg = [0.5 0.5 scale];   % position of Fig 1
set(figure(1),'position',[sw0*pg(1) sh0*pg(2)  sw0*pg(3) sh0*pg(4)])

% Figure 2
scale = [1 .83]*0.5;     % horiz/vert size
pg = [0.5 0.0 scale];	% position of Fig 2
set(figure(2),'position',[sw0*pg(1) sh0*pg(2)  sw0*pg(3) sh0*pg(4)])

%%% CONTROL PANEL = figure(500)
scale = [.31 1];     % horiz/vert size
pg = [0.18 0.05 scale];	% position of GUI window
guiCtrl = figure(500);clf
set(guiCtrl,'Resize','on','Units','pixels','Position',[sw0*pg(1) sh0*pg(2)  sw0*pg(3) sh0*pg(4)],'Visible','off',...
    'MenuBar','none','name','GUI_mcxyz','NumberTitle','off','UserData',0);

%% %% MESSAGES
defaultBackground = get(0,'defaultUicontrolBackgroundColor');
set(guiCtrl,'Color',defaultBackground);
set(guiCtrl,'Visible','on');

% MESSAGE box
note1 = 'Ready to make a tissue.';
infoLabel = uicontrol('Parent',guiCtrl,'Style','text','String',note1,'FontUnits',...
    'normalized','FontSize',.12,'Units','normalized','Position',[0 .00 1.0 .2],'BackgroundColor',[.7 1 1]);
set(infoLabel,'String',note1);

%%%%%% 
% CONTROL BUTTONS
%%%

% MAKETISSUE  
xB = 0; yB=0.72;     % position of button
wB = .45; hB = .05;   % size of button
name = 'make tissue';
callMAKETISSUE = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',name,'BackgroundColor',clr(1,:),...
    'FontUnits','normalized','FontSize',.3,'fontweight','bold','Units','normalized','Position',[xB yB wB hB],...
    'Callback',{@MAKETISSUE}); 

% GOMCXYZ to run mcxyz
xB = 0; yB=0.65;     % position of button
GOcolor = 'g'; 
name = 'GO mcxyz';
callGOMCXYZ = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',name,'BackgroundColor',clr(3,:),...
    'FontUnits','normalized','FontSize',.3,'fontweight','bold','Units','normalized','Position',[xB yB wB hB],...
    'Callback',{@GOMCXYZ}); 

% GOMCXYZBKGD to run mcxyz
xB = 0; yB=0.58;     % position of button
GOcolor = 'g'; 
name = 'GO mcxyz BKGD';
callGOMCXYZBKGD = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',name,'BackgroundColor',clr(3,:),...
    'FontUnits','normalized','FontSize',.3,'fontweight','bold','Units','normalized','Position',[xB yB wB hB],...
    'Callback',{@GOMCXYZBKGD}); 

% SHOWMCXYZ
xB = 0; yB=0.51;     % position of button
SHOWcolor = [.7 .7 1]; 
name = 'SHOW mcxyz';
callSHOWMCXYZ = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',name,'BackgroundColor',clr(3,:),...
    'FontUnits','normalized','FontSize',.3,'fontweight','bold','Units','normalized','Position',[xB yB wB hB],...
    'Callback',{@SHOWMCXYZ}); 

% reload BUTTON
xB = 0; yB=0.92;     % position of button
name = 'reload simulation';
callRELOAD = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',name,'BackgroundColor',clr(1,:),...
    'FontUnits','normalized','FontSize',.3,'fontweight','bold','Units','normalized','Position',[xB yB wB hB],...
    'Callback',{@CALLreload}); 

% test BUTTON
xB = 0; yB=0.25;     % position of button
THERMALcolor = [1 .7 .7]; 
name = 'test';
callTEST = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',name,'BackgroundColor',clr(3,:),...
    'FontUnits','normalized','FontSize',.3,'fontweight','bold','Units','normalized','Position',[xB yB wB hB],...
    'Callback',{@CALLtest}); 

% BUTTON to reset gui
xB = 0.82; yB=0.01;     % position of button
wB = .15; hB = .04;   % size of button
name = 'Reset';
RESET = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',name,'BackgroundColor','w',...
    'FontUnits','normalized','FontSize',.50,'Fontweight','bold','Units','normalized',...
    'Position',[xB yB wB hB],'Callback',{@callRESET});

% QUIT BUTTON
xB = 0.92; yB=0.96;     % position of button
wB = .06; hB = .03;   % size of button
name = 'X';
callQUIT = uicontrol('Parent',guiCtrl,'Style','pushbutton','String',name,'BackgroundColor',[1 1 1]*0.8,...
    'FontUnits','normalized','FontSize',.4,'fontweight','bold','Units','normalized',...
    'Position',[xB yB wB hB],'Callback',{@CALLquit}); 

%%%%%% 
% Editable text and table
%%%
xB = 0; xB2 = 0.13; yB=0.905; dyB = .04; % position of labels & textboxes 
wB = .13; wB2 = .4;hB = .03;  % size of button
fz = .4;

% MC controls 
calltextbox_MCcontrols = uicontrol('Parent',guiCtrl,'Style','text','String',['Monte Carlo controls'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[0.55 0.89  0.4 hB]);

s = 'zfocus is relative to zs';
calltextbox_msg1 = uicontrol('Parent',guiCtrl,'Style','text','String',s,...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[0.55 0.42  0.4 hB]);

% pick a name
j=1;
calltextbox_name = uicontrol('Parent',guiCtrl,'Style','text','String',['simname'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
entersimname = uicontrol('Parent',guiCtrl,'Style','edit','String',simname,'BackgroundColor',[1 .5 .5],...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB2 hB],...
    'fontsize',12,'Callback',{@getsimname});

% pick a tissue
j=2;
calltextbox_tissue = uicontrol('Parent',guiCtrl,'Style','text','String',['tissuename'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
entertissuename = uicontrol('Parent',guiCtrl,'Style','edit','String',tissuename,'BackgroundColor',[1 1 1],...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB2 hB],...
    'fontsize',12,'Callback',{@gettissuename});

% pick a wavelength
j=3;
calltextbox_name = uicontrol('Parent',guiCtrl,'Style','text','String',['wavelength'],...
    'FontUnits','normalized','FontSize',fz,'Units','normalized','Position',[xB yB-dyB*j wB hB]);
enterwavelength = uicontrol('Parent',guiCtrl,'Style','edit','String','wavelength [nm]','BackgroundColor',[1 .5 .5],...
    'Min',0,'Max',1,'UserData',[],'Units','normalized','Position',[xB2 yB-dyB*j wB2 hB],...
    'fontsize',12,'Callback',{@getwavelength});

% editable TABLE of Monte Carlo parameters
cnames = 'parameters';
Nmc = length(mc); % # of variables in table (some are empty)
u = java_array('java.lang.String',Nmc);
dat = zeros(Nmc,1);
for i=1:Nmc
    u(i) = java.lang.String(mc(i).name);
    dat(i) = mc(i).param;
end
rnames = cell(u);
pos = [.54 .45 .45 .40];
t = uitable('Parent',guiCtrl,'Data',dat,'ColumnName',cnames,... 
            'RowName',rnames,'Units','normalized','Position',pos,...
            'ColumnEditable',true(1,20));


%%%%%%%%%%%%%%
% CALLBACKS
%%%

%--------------------------------------------------------------------------
    function getsimname(entersimname,eventdata)
        simname = get(entersimname,'string');
        fprintf('name = %s\n',simname);
        set(entersimname,'Background',[.7 1 .7])
        pause(0.5)
        set(entersimname,'Background','w')
        set(callMAKETISSUE,'backgroundColor',clr(2,:)) 
        set(callGOMCXYZ,'backgroundColor',clr(1,:)) 
        set(callGOMCXYZBKGD,'backgroundColor',clr(1,:)) 
        set(callSHOWMCXYZ,'backgroundColor',clr(1,:)) 
    end

%--------------------------------------------------------------------------
    function gettissuename(entertissuename,eventdata)
        tissuename = get(entertissuename,'string');
        fprintf('name = %s\n',tissuename);
        set(entertissuename,'Background',[.7 1 .7])
        pause(0.5)
        set(entertissuename,'Background','w')
        set(callMAKETISSUE,'backgroundColor',clr(2,:)) 
        set(callGOMCXYZ,'backgroundColor',clr(1,:)) 
        set(callGOMCXYZBKGD,'backgroundColor',clr(1,:)) 
        set(callSHOWMCXYZ,'backgroundColor',clr(1,:)) 
    end

%--------------------------------------------------------------------------
    function getwavelength(enterwavelength,eventdata)
        u = round(get(enterwavelength,'string'));
        ii = find(u>=48 & u<=57);
        nm = str2num(char(u(ii)));
        fprintf('wavelength = %d nm\n',nm);
        set(enterwavelength,'Background',[.7 1 .7])
        pause(0.5)
        set(enterwavelength,'Background','w','string',[num2str(nm) ' nm'])
        set(callMAKETISSUE,'backgroundColor',clr(2,:)) 
        set(callGOMCXYZ,'backgroundColor',clr(1,:)) 
        set(callGOMCXYZBKGD,'backgroundColor',clr(1,:)) 
        set(callSHOWMCXYZ,'backgroundColor',clr(1,:)) 
    end

%--------------------------------------------------------------------------
    function MAKETISSUE(callMAKETISSUE,eventdata)
        if strcmp(simname,'pick_a_name')
            s = sprintf('First, pick a name\nfor the tissue created.\n(always hit return to confirm entry)');
            disp(s)
        	set(infoLabel,'String',s);
        else
            set(callMAKETISSUE,'backgroundColor',clr(2,:)) 
            pause(0.5)
            % read table t --> mc.
            u = get(t,'Data');
            for i=1:Nmc
                mc(i).param = u(i);
            end
            % save params
            mkdir(sprintf('../simsLibrary/%s',simname))
            fid = fopen(sprintf('../simsLibrary/%s/params.dat',simname),'w');
            for i=1:Nmc
                fprintf(fid,'%0.5f\n',u(i))
            end
            fclose(fid);
            %
            set(callMAKETISSUE,'backgroundColor',clr(1,:)) 
            set(callGOMCXYZ,'backgroundColor',clr(1,:)) 
            set(callGOMCXYZBKGD,'backgroundColor',clr(1,:)) 
            set(callSHOWMCXYZ,'backgroundColor',clr(1,:)) 
            pause(0.2)
            fprintf('\n\n\n\nmake tissue @ %d nm\n',nm)
            disp('go maketissue');
            maketissue(mc,tissuename,simname,nm);
            figure(guiCtrl)
            set(infoLabel,'String',...
                sprintf('tissue "%s" created.\nReady for GO mcxyz.',simname)); 
            set(callMAKETISSUE,'backgroundColor','w')
            set(callGOMCXYZ,'backgroundColor',clr(2,:))
            set(callGOMCXYZBKGD,'backgroundColor',clr(2,:))
            GOflag=1;
        end
    end

%--------------------------------------------------------------------------
    function GOMCXYZ(callGOMCXYZ,eventdata)
        if GOflag
            set(callGOMCXYZ,'backgroundColor',clr(4,:))        
            pause(.1)
            if ismac
                % run mcxyz on MAC
                fprintf('\nrunning mcxyz...\n')
                set(infoLabel,'String','running mcxyz');
                pause(0.1)
                cd(homedir)
                cd(sprintf('../simsLibrary/%s',simname))
                cmd = sprintf('%s/gomcxyz %s',homedir,simname);
                system(cmd);
                cd(homedir)
            else
                % run mcxyz on WINDOWS
                fprintf('\nrunning mcxyz...\n')
                set(infoLabel,'String','running mcxyz');
                pause(0.1)
                cd(homedir)
                cd(sprintf('../simsLibrary/%s',simname))
                cmd = sprintf('%s/mcxyz_W7.exe %s',homedir,simname);
                system(cmd);
                cd(homedir)
            end
            set(callGOMCXYZ,'backgroundColor',clr(1,:))
            set(infoLabel,'String','mcxyz done.');
            GOflag==0;
            SHOWflag = 1;
            set(callSHOWMCXYZ,'backgroundColor',clr(2,:))
        else
            fprintf('Not ready for GO mcxyz.\nUse "make tissue"')
            set(infoLabel,'String',...
                sprintf('Not ready for GO mcxyz.\nUse "make tissue"'));
        end
    end

%--------------------------------------------------------------------------
    function GOMCXYZBKGD(callGOMCXYZBKGD,eventdata)
        if GOflag==1
            set(callGOMCXYZBKGD,'backgroundColor',clr(4,:))        
            pause(.1)
            if ismac
                % run mcxyz on maxOSX
                    fprintf('\nrunning mcxyz in the background...\n')
                    set(infoLabel,'String','running mcxyz in background');
                    pause(0.1)
                    cd(homedir)
                    cd(sprintf('../simsLibrary/%s',simname))
                    cmd = sprintf('%s/gomcxyz %s &',homedir,simname);
                    system(cmd);
                    !ps
                    cd(homedir)
             else
                % run mcxyz on WINDOWS
                    fprintf('\nrunning mcxyz in the background...\n')
                    set(infoLabel,'String','running mcxyz in background');
                    pause(0.1)
                    cd(homedir)
                    cd(sprintf('../simsLibrary/%s',simname))
                    cmd = sprintf('%s/mcxyz_W7.exe %s &',homedir,simname);
                    system(cmd);
                    cd(homedir)
            end
            set(callGOMCXYZBKGD,'backgroundColor',clr(1,:))
        	fprintf('mcxyz done.\n\n')
            set(infoLabel,'String','mcxyz done.');
            GOflag==0;
            SHOWflag = 1;
            set(callSHOWMCXYZ,'backgroundColor',clr(2,:))
        else
            set(infoLabel,'String',...
                sprintf('Not ready for GO mcxyz.\nUse "make tissue"'));
        end
    end

%--------------------------------------------------------------------------
    function SHOWMCXYZ(callSHOWMCXYZ,eventdata)
        fprintf('\nshow mcxyz output\n')
        %try
            tissue = makeTissueList(nm,1); % 1=PRINTON
            [Qyxz Fyxz] = showmcxyz(simname,nm,mc,tissue);
            set(callSHOWMCXYZ,'backgroundColor',clr(1,:))
            s = 'fluence and deposition';
            set(infoLabel,'String',s);
%         catch
%             fprintf('sorry, %s not ready\n',simname)
%         end
        figure(guiCtrl)
        sfig1 = 'A';
    end


%--------------------------------------------------------------------------
    function CALLreload(callRELOAD,eventdata)
        fprintf('\nRELOAD: \n')
            set(callRELOAD,'backgroundColor',clr(2,:))
            pause(.5)
            commandwindow
            
            D = dir('../simsLibrary/');
            for i=1:length(D)
                if ~strcmp(D(i).name(1),'.')
                    fprintf('%d\t%s\n',i,D(i).name)
                end
            end
            
            try
                j = input('pick a simulation by # ');
                simname = D(j).name;
            catch
                disp('no change')
                return
            end            
            set(entersimname,'string',simname)
            
          	load(sprintf('../simsLibrary/%s/tissuename.mat',simname))
                % --> tissue().xx, nm
            set(enterwavelength,'string',[num2str(nm) ' nm'])
            set(entertissuename,'string',tissuename)
            
            u = load(sprintf('../simsLibrary/%s/params.dat',simname));
            Nu = length(u);
            for i=1:Nu
                dat(i) = u(i);                
            end            
            for i=1:Nu
                mc(i).param = u(i);
            end
            %pos = [.54 .45 .45 .45];
            t = uitable('Parent',guiCtrl,'Data',dat,'ColumnName',cnames,... 
                        'RowName',rnames,'Units','normalized','Position',pos,'ColumnEditable',...
                        true(1,20));
            GOflag = 1;
            
        set(callRELOAD,'backgroundColor',clr(1,:))
        set(entersimname,'backgroundColor',clr(1,:))
        set(entertissuename,'backgroundColor',clr(1,:))
        set(enterwavelength,'backgroundColor',clr(1,:))
        set(callGOMCXYZ,'backgroundColor',clr(2,:))
        set(callGOMCXYZBKGD,'backgroundColor',clr(2,:))
        set(callSHOWMCXYZ,'backgroundColor',clr(2,:))
        figure(guiCtrl)
    end

%--------------------------------------------------------------------------
    function CALLtest(callTEST,eventdata)
        fprintf('\nTESTING: \n')
            set(callTEST,'backgroundColor',clr(2,:))
            pause(.5)
            disp('TEST CODE:')
            %%%
            % place test code here:
            %disp('no test code')
          
            CH = 1;
            switch CH
                case 1
                    !ps
                    nm
                case 2
                    fprintf('simulation name = %s\n',simname)
                    ls(sprintf('../simsLibrary/%s',simname))
                    fprintf('GOflag = %d\n',GOflag)

                    zfocus = mc(16).param
                    radius = mc(17).param
                    ang=atan(zfocus/radius)
                    deg = ang*180/pi
                    NA = sin(ang)
                    angi = ang*1.4
                    degi = angi*180/pi
                    NAi = sin(angi)

                    zfocus
                case 3
                    fprintf('wavelength = %0.0f nm\n',nm)
            end % switch
            
            % end code.
            %%%
            
        set(callTEST,'backgroundColor',clr(1,:))
        figure(guiCtrl)
    end

%--------------------------------------------------------------------------
% callback function for QUIT button
    function CALLquit(callQUIT,eventdata)
        disp('Quit.')
        close(1)
        close(2)
        close(500)
    end

%--------------------------------------------------------------------------
% callback function for RESET
    function callRESET(RESET,eventdata)
        close(500)
        eval(gui)
    end


    disp('ready')
end % ends optogensim function
