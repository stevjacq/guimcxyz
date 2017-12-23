function reportHmci(simname)
% function reportHmci(simname)
%   Lists the values of the input file simname_H.mci.
%   Updated Feb 8, 2017. slj, adding boundaryflag B(10) (see s(10).s)

home
fid = fopen(sprintf('../simsLibrary/%s/%s_H.mci',simname,simname),'r');
B = fscanf(fid,'%f');
fclose(fid);

s(1).s = 'time_min';
s(2).s = 'Nx';
s(3).s = 'Ny';
s(4).s = 'Nz';
s(5).s = 'dx';
s(6).s = 'dy';
s(7).s = 'dz';
s(8).s = 'mcflag';
s(9).s = 'boundary';
s(10).s = 'xs';
s(11).s = 'ys';
s(12).s = 'zs';
s(13).s = 'xfocus';
s(14).s = 'yfocus';
s(15).s = 'zfocus';
s(16).s = 'radius';
s(17).s = 'waist';
s(18).s = 'zsurf';
s(19).s = 'launch';
s(20).s = 'ux0';
s(21).s = 'uy0';
s(22).s = 'uz0';
s(23).s = 'Nt';

for i=1:23
    disp(sprintf('%d\t%10s = %0.4f',i,s(i).s,B(i)))
end

if 0 % 1 = show optical properties of tissues
for j=1:B(23)
    i=i+1;
    disp(sprintf('---'))
    disp(sprintf('%d\tmua = %0.4f',i,B(i)))
    i=i+1;
    disp(sprintf('%d\tmus = %0.4f',i,B(i)))
    i=i+1;
    disp(sprintf('%d\tg   = %0.4f',i,B(i)))
end
end
