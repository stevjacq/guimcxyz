function [Qyxz F] = showmcxyz(simname,nm,mc,tissue)
%   Looks at simname_F.bin, created by mcxyz.c 
%   where simname is the name of the run: simname_T.bin, simname_H.mci
%   Makes figures:
%       simname_tissue.jpg   = tissue structure (shows tissue types)
%       simname_Fzx.jpg      = fluence rate vs z,x
%       simname_Fzy.jpg      = fluence rate vs z,y
%   Uses:
%       simname_H.mci    = input file from maketissue.m
%       simname_T.bin    = tissue input file from maketissue.m
%       simname_F.bin    = fluence rate output from Monte Carlo
%       reportH_mci.m   = lists input parameters in simname_H.mci
%       makecmap.m      = makes colormap for tissue types
%       makec2f.m       = makes colormap for fluence rate
%
%   This example sets simname = 'skinvessel'.
%
% 7/feb/2017, add boundaryflag (see A(10)).
% 1/june/2017 , no major changes, just clean up display outputs.
% Steven L Jacques
home
format compact
commandwindow

SAVEPICSON = 0;
if SAVEPICSON
    sz = 10; fz = 7; fz2 = 5; % to use savepic.m
else
    sz = 14; fz = 12; fz2 = 12; % for screen display
end


%%%% USER CHOICES <---------- you must specify -----


disp(sprintf('------ mcxyz %s -------',simname))

% Load header file
filename = sprintf('../simsLibrary/%s/%s_H.mci',simname,simname);
disp(['loading ' filename])
fid = fopen(filename, 'r');
A = fscanf(fid,'%f',[1 Inf])';
fclose(fid);

%% parameters
time_min = A(1);
Nx = A(2);
Ny = A(3);
Nz = A(4);
dx = A(5);
dy = A(6);
dz = A(7);
mcflag = A(8);
boundaryflag = A(9);
xs = A(10);
ys = A(11);
zs = A(12);
xfocus = A(13);
yfocus = A(14);
zfocus = A(15);
radius = A(16);
waist = A(17);
zsurf = A(18);
launchflag = A(19);
ux0 = A(20);
uy0 = A(21);
uz0 = A(22);
Nt = A(23);
j = 23;
for i=1:Nt
    j=j+1;
    muav(i,1) = A(j);
    j=j+1;
    musv(i,1) = A(j);
    j=j+1;
    gv(i,1) = A(j);
end

reportHmci(simname)

%% Load Fluence rate F(y,x,z) 
filename = sprintf('../simsLibrary/%s/%s_F.bin',simname,simname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'float');
    fclose(fid);
toc
F = reshape(Data,Ny,Nx,Nz); % F(y,x,z)

%%
% Load tissue structure in voxels, T(y,x,z) 
filename = sprintf('../simsLibrary/%s/%s_T.bin',simname,simname);
disp(['loading ' filename])
tic
    fid = fopen(filename, 'rb');
    [Data count] = fread(fid, Ny*Nx*Nz, 'uint8');
    fclose(fid);
toc
T = reshape(Data,Ny,Nx,Nz); % T(y,x,z)
clear Data

%% Rd
Rd = load(sprintf('../simsLibrary/%s/%s_Rd.dat',simname,simname))

%%
x = ([1:Nx]-Nx/2-1/2)*dx;
y = ([1:Ny]-Ny/2-1/2)*dx;
z = ([1:Nz]-1/2)*dz;
ux = [2:Nx-1];
uy = [2:Ny-1];
uz = [2:Nz-1];
zmin = min(z);
zmax = max(z);
zdiff = zmax-zmin;
xmin = min(x);
xmax = max(x);
xdiff = xmax-xmin;

%% Look at structure, Tzx
Tzx = reshape(T(Ny/2,:,:),Nx,Nz)'; 
tissue = makeTissueList(nm,1); % 0 = PRINTON
Nt = length(tissue);

figure(1);clf
imagesc(x(ux),z(uz),Tzx(uz,ux),[1 Nt])
hold on
cmap = makecmap(Nt);
colormap(cmap)
colorbar
set(gca,'fontsize',sz)
set(colorbar,'fontsize',1)
xlabel('x [cm]')
ylabel('z [cm]')
title('Tissue','fontweight','normal','fontsize',fz)
for i=1:Nt
    yy = zmin + (Nt-i)/(Nt-1)*zdiff;
    text(xmax*1.2,yy, sprintf('%d %s',i,tissue(i).name),...
        'fontsize',fz2)
end

% draw launch
N = 20; % # of beam rays drawn
clrsym = [1 .3 .3];%'c--';
switch mcflag
    case 0 % uniform
        for i=0:N
            if zfocus==inf
                plot([xs+(-radius + 2*radius*i/N) xfocus],[zs max(z)],'r-')
            else
                plot([xs+(-radius + 2*radius*i/N) xfocus],zs+[0 zfocus],'color',clrsym)
            end
        end

    case 1 % Gaussian
        for i=0:N
            plot([xs+(-radius + 2*radius*i/N) xfocus],zs+[0 zfocus],'color',clrsym)
        end

    case 2 % iso-point
        for i=1:N
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'color',clrsym)
        end
        
    case 3 % rectangle
        zz = max(z);
        for i=0:N
            xx = -radius + 2*radius*i/N;
            plot([xx xx],[zs zz],'color',clrsym)
        end
end
axis equal image
axis([min(x) max(x) min(z) max(z)])

if SAVEPICSON
    name = sprintf('%s_tissue.png',simname);
    savepic(1,[4 3],name)
end

%% Qyxz(y,x,z) [J/cm3 per J delivered] = [1/cm^3]
for i=1:length(tissue)
    muav(i) = tissue(i).mua;
end
Qyxz = reshape(F.*muav(T),Ny,Nx,Nz);




%% Look at Fluence Fzx @ launch point
Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source

figure(2);clf
subplot(1,2,1)
imagesc(x,z,log10(Fzx))
hold on
text(max(x),min(z)-0.05*max(z),'log_{10}( \phi )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('Fluence \phi [J/cm^2/J.delivered] ','fontweight','normal','fontsize',fz)
colormap(makec2f)
axis equal image
text(x(1)-x(end)*.4,0,'A','fontsize',sz)
%axis([min(x) max(x) min(z) max(z)])
text(min(x)-0.2*max(x),min(z)-0.15*max(z),sprintf('runtime = %0.1f min',time_min),...
    'fontsize',fz2)

if SAVEPICSON
    name = sprintf('%s_Fzx.png',simname);
    savepic(10,[4 3],name)
end


%% look Azx
Fzx = reshape(F(Ny/2,:,:),Nx,Nz)'; % in z,x plane through source
mua = muav(reshape(T(Ny/2,:,:),Nx,Nz)');
Azx = Fzx.*mua;

% zbulb = tiss(1).param-0.0400; % through center of bulb
% izbulb = round(zbulb/dz);
% Tbx = mean(Azx(izbulb+[-2:2],:))/4.2; % avg 5 zbins @ zbulb 
% Tez = mean(Azx(:,30:40),2)/4.2; % average over range of x's

W = 0.65;
rhoCp = 4.2*W;

figure(2)
subplot(1,2,2)
imagesc(x,z,log10(Azx/rhoCp))
hold on
%plot([x(1) x(end)],[1 1]*log10(Azx(izbulb,:)),'w--')
text(max(x),min(z)-0.05*max(z),'log_{10}( \DeltaT )','fontsize',fz)
colorbar
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
title('\DeltaT [\circC/J.delivered] ','fontweight','normal','fontsize',fz)
text(0,-max(z)/6,sprintf('water content = %0.2f',W))
colormap(makec2f)
axis equal image
text(x(1)-x(end)*.4,0,'B','fontsize',sz)

    
    %axis([min(x) max(x) min(z) max(z)])
% text(min(x)-0.2*max(x),min(z)-0.08*max(z),sprintf('runtime = %0.1f min',time_min),...
%     'fontsize',fz2)

% ii=find(z>=0.275); iz1 = ii(1);
% figure(4)
% plot([-.2 .2],[1 1]*z(iz1),'--','color',[1 1 1]*0.99)

if SAVEPICSON
    name = sprintf('%s_Azx.png',simname);
    savepic(2,[4 3],name)
end

drawnow


