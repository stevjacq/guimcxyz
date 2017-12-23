function maketissue(mc,tissuename,simname,nm)
% maketissue_18apr17.m
% maketissue.m
%   Creates a cube of optical property pointers,T(y,x,z), saved in
%       simname_T.bin = a tissue structure file
%   which specifies a complex tissue for use by mcxyz.c.
%
%   Also prepares a listing of the optical properties at chosen wavelength
%   for use by mcxyz.c, [mua, mus, g], for each tissue type specified
%   in simname_T.bin. This listing is saved in
%       simname_H.mci = the input file for use by mcxyz.c.
%
%   Will generate a figure illustrating the tissue with its various
%   tissue types and the beam being launched.
%
%   Uses
%       makeTissueList.m
%
%   To use, 
%       1. Prepare makeTissueList.m so that it contains the tissue
%   types desired.
%       2. Specify the USER CHOICES.
%       2. Run this program, maketissue.m.
%
%   Note: mcxyz.c can use optical properties in cm^-1 or mm^-1 or m^-1,
%       if the bin size (binsize) is specified in cm or mm or m,
%       respectively.
%
%  Steven L. Jacques. updated Aug 21, 2014.
%       

format compact
homedir = '~/Documents/GitHub/mcxyzGUI';
path(path,homedir)
cd(homedir)

dirname     = sprintf('~/Documents/GitHub/simsLibrary/%s',simname);
mkdir(dirname)


%%% USER CHOICES %%%%%%%% <-------- You must set these parameters ------
SAVEON      = 1;        % 1 = save simname_T.bin, simname_H.mci 
                        % 0 = don't save. Just check the program.

i=1;                        
time_min	= mc(i).param;
i=i+1;
Nx          = mc(i).param;
i=i+1;
Ny          = mc(i).param;
i=i+1;
Nz          = mc(i).param;
i=i+1;
dx          = mc(i).param;
i=i+1;
dy          = mc(i).param;
i=i+1;
dz          = mc(i).param;
i=i+1;
mcflag      = mc(i).param;
i=i+1;
boundaryflag= mc(i).param;
i=i+1;
xs          = mc(i).param;
i=i+1;
ys          = mc(i).param;
i=i+1;
zs          = mc(i).param;
i=i+1;
xfocus      = mc(i).param;
i=i+1;
yfocus      = mc(i).param;
i=i+1;
zfocus      = mc(i).param;
i=i+1;
radius      = mc(i).param;
i=i+1;
waist       = mc(i).param;
i=i+1;
zsurf       = mc(i).param;
i=i+1;
launchflag  = mc(i).param;
i=i+1;
ux0         = mc(i).param;
i=i+1;
uy0         = mc(22).param;
i=i+1;
uz0         = mc(i).param;

% Monte Carlo flags:
% mcflag launch: 0 = uniform beam, 1 = Gaussian, 2 = isotropic pt. 
%   3 = rectangular beam (use xfocus,yfocus for x,y halfwidths)
% launchflag 0 = let mcxyz.c calculate launch trajectory
%   1 = manually set launch vector.
% boundaryflag 0 = no boundaries, 1 = escape at boundaries
%   2 = escape at surface only. No x, y, bottom z boundaries



%%%%%%%%%% 
% Prepare Monte Carlo 
%%%

% Create tissue properties
PRINTON = 1
tissue = makeTissueList(nm,PRINTON); % also --> global tissue(1:Nt).s
Nt = length(tissue);
for i=1:Nt
    muav(i)  = tissue(i).mua;
    musv(i)  = tissue(i).mus;
    gv(i)    = tissue(i).g;
end

% Specify Monte Carlo parameters    
x  = ([1:Nx]'-Nx/2)*dx;
y  = ([1:Ny]'-Ny/2)*dy;
z  = [1:Nz]'*dz;
zmin = min(z);
zmax = max(z);
xmin = min(x);
xmax = max(x);

if isinf(zfocus), zfocus = 1e12; end

%%%%%%
% CREATE TISSUE STRUCTURE T(y,x,z)
%   Create T(y,x,z) by specifying a tissue type (an integer)
%   for each voxel in T.
%
%   Note: one need not use every tissue type in the tissue list.
%   The tissue list is a library of possible tissue types.

T = getT(tissuename,mc);
disp(tissuename)

%%
if SAVEON
    
    filename = sprintf('../simsLibrary/%s/tissuename.mat',simname)
    save(filename,'tissuename','nm');
    
    tic
    % convert T to linear array of integer values, v(i)i = 0;
    v = uint8(reshape(T,Ny*Nx*Nz,1));

    %% WRITE FILES
    % Write simname_H.mci file
    %   which contains the Monte Carlo simulation parameters
    %   and specifies the tissue optical properties for each tissue type.
    commandwindow
    disp(sprintf('--------create %s --------',simname))
    filename = sprintf('../simsLibrary/%s/%s_H.mci',simname,simname);
    fid = fopen(filename,'w');
        % run parameters
        fprintf(fid,'%0.2f\n',time_min);
        fprintf(fid,'%d\n'   ,Nx);
        fprintf(fid,'%d\n'   ,Ny);
        fprintf(fid,'%d\n'   ,Nz);
        fprintf(fid,'%0.4f\n',dx);
        fprintf(fid,'%0.4f\n',dy);
        fprintf(fid,'%0.4f\n',dz);
        % launch parameters
        fprintf(fid,'%d\n'   ,mcflag);
        fprintf(fid,'%d\n'   ,boundaryflag);
        fprintf(fid,'%0.4f\n',xs);
        fprintf(fid,'%0.4f\n',ys);
        fprintf(fid,'%0.4f\n',zs);
        fprintf(fid,'%0.4f\n',xfocus);
        fprintf(fid,'%0.4f\n',yfocus);
        fprintf(fid,'%0.4f\n',zfocus);
        fprintf(fid,'%0.4f\n',radius);
        fprintf(fid,'%0.4f\n',waist);
        fprintf(fid,'%0.4f\n',zsurf);
        fprintf(fid,'%d\n'   ,launchflag); % if manually setting ux,uy,uz
        fprintf(fid,'%0.4f\n',ux0); 
        fprintf(fid,'%0.4f\n',uy0);
        fprintf(fid,'%0.4f\n',uz0);
        % tissue optical properties
        fprintf(fid,'%d\n',Nt);
        for i=1:Nt
            fprintf(fid,'%0.4f\n',muav(i));
            fprintf(fid,'%0.4f\n',musv(i));
            fprintf(fid,'%0.4f\n',gv(i));
        end
    fclose(fid);

    reportHmci(simname)

    %% write simname_T.bin file
    filename = sprintf('../simsLibrary/%s/%s_T.bin',simname,simname);
    disp(['create ' filename])
    fid = fopen(filename,'wb');
    fwrite(fid,v,'uint8');
    fclose(fid);

    toc
end % SAVEON


%% Look at structure of Tzx at iy=Ny/2
Txzy = shiftdim(T,1);   % Tyxz --> Txzy
Tzx  = Txzy(:,:,Ny/2)'; % Tzx

%%
figure(1); clf
sz = 12;  fz = 10; 
imagesc(x,z,Tzx,[1 Nt])
hold on
set(gca,'fontsize',sz)
xlabel('x [cm]')
ylabel('z [cm]')
colorbar
cmap = makecmap(Nt);
colormap(cmap)
set(colorbar,'fontsize',1)
% label colorbar
zdiff = zmax-zmin;
%%%

for i=1:Nt
    yy = (Nt-i)/(Nt-1)*Nz*dz;
    text(max(x)*1.2,yy, tissue(i).name,'fontsize',fz)
end

text(xmax,zmin - zdiff*0.06, 'Tissue types','fontsize',fz)
axis equal image
axis([xmin xmax zmin zmax])

%%% draw launch
N = 20; % # of beam rays drawn
switch mcflag
    case 0 % uniform
        for i=0:N
            if zfocus==inf
                plot([xs+(-radius + 2*radius*i/N) xfocus],[zs 1e3],'r-')
            else
                plot([xs+(-radius + 2*radius*i/N) xfocus],zs+[0 zfocus],'r-')
            end
        end

    case 1 % Gaussian
        for i=0:N
            plot([xs+(-radius + 2*radius*i/N) xfocus],zs+[0 zfocus],'r-')
        end

    case 2 % iso-point
        for i=1:N
            th = (i-1)/19*2*pi;
            xx = Nx/2*cos(th) + xs;
            zz = Nx/2*sin(th) + zs;
            plot([xs xx],[zs zz],'r-')
        end
        
    case 3 % rectangle
        zz = max(z);
        for i=1:N
            xx = -radius + 2*radius*i/20;
            plot([xx xx],[zs zz],'r-')
        end
end

disp('done')

