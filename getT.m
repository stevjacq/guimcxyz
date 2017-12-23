function T = getT(tissuename,mc)
% function T = getT(tissuename,mc)
%   getT; simply lists the available tissue models.
% The available tissue models are defined in this routine.
% You can add tissue models.

% Just list tissue types
if nargin==0
    tiss(1).s = 'skin_vessel';
    tiss(2).s = 'skin_model';
    tiss(3).s = 'soft_tissue';
    for i=1:length(tiss)
        fprintf('%d\t%s\n',i,tiss(i).s)
    end
    T = 0;
    return;
end

% Create T(y,x,z)
Nx          = mc(2).param;
Ny          = mc(3).param;
Nz          = mc(4).param;
dx          = mc(5).param;
dy          = mc(6).param;
dz          = mc(7).param;
zsurf       = mc(18).param;

x  = ([1:Nx]'-Nx/2)*dx;
y  = ([1:Ny]'-Ny/2)*dy;
z  = [1:Nz]'*dz;

T = double(zeros(Ny,Nx,Nz)); 

switch tissuename
    case 'skin_vessel'
        T = T + 4;      % fill background with skin (dermis)
        for iz=1:Nz % for every depth z(iz)
            % air
            if iz<=round(zsurf/dz)
                T(:,:,iz) = 2; 
            end
            % epidermis (60 um thick)
            if iz>round(zsurf/dz) & iz<=round((zsurf+0.0060)/dz)
                T(:,:,iz) = 5; 
            end
            % blood vessel @ xc, zc, radius, oriented along y axis
            xc      = 0;            % [cm], center of blood vessel
            zc      = Nz/2*dz;     	% [cm], center of blood vessel
            vesselradius  = 0.0100;      	% blood vessel radius [cm]
            for ix=1:Nx
                    xd = x(ix) - xc;	% vessel, x distance from vessel center
                    zd = z(iz) - zc;   	% vessel, z distance from vessel center                
                    r  = sqrt(xd^2 + zd^2);	% r from vessel center
                    if (r<=vesselradius)    % if r is within vessel
                        T(:,ix,iz) = 3; % blood
                    end

            end %ix
        end % iz

    case 'skin_model'
        T = T + 4;
        iz1 = round(zsurf/dz);
        iz2 = round((zsurf+0.0050)/dz);
        iz3 = round((zsurf+0.0060)/dz);
        ibc = [iz1 iz2 iz3]
        for iz=1:iz1 % air            
                T(:,:,iz) = 2; 
        end
        for iz=iz1+1:iz2 % s.epidermis            
                T(:,:,iz) = 5; 
        end
        for iz=iz3 % basal           
                T(:,:,iz) = 9; 
        end
        
    case 'soft_tissue'
        T = T + 10;
        iz1 = round(zsurf/dz);
        for iz=1:iz1
            T(:,:,iz) = 1;
        end
        
end % switch