function reportHmci(myname)
% function reportHmci(myname)
%   Lists the values of the input file  myname_H.mci.
%   Updated Feb 8, 2017. slj, adding boundaryflag B(10) (see s(10).s)

home
fid = fopen(sprintf('%s_H.mci',myname),'r');
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
s(9).s = 'launch';
s(10).s = 'boundary';
s(11).s = 'xs';
s(12).s = 'ys';
s(13).s = 'zs';
s(14).s = 'xfocus';
s(15).s = 'yfocus';
s(16).s = 'zfocus';
s(17).s = 'ux0';
s(18).s = 'uy0';
s(19).s = 'uz0';
s(20).s = 'radius';
s(21).s = 'waist';
s(22).s = 'Nt';

for i=1:22
    disp(sprintf('%d\t%10s = %0.4f',i,s(i).s,B(i)))
end

for j=1:B(22)
    i=i+1;
    disp(sprintf('---'))
    disp(sprintf('%d\tmua = %0.4f',i,B(i)))
    i=i+1;
    disp(sprintf('%d\tmus = %0.4f',i,B(i)))
    i=i+1;
    disp(sprintf('%d\tg   = %0.4f',i,B(i)))
end
