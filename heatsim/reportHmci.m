function reportHmci(directoryPath,myname)
% function reportHmci(myname)
%   Lists the values of the input file  myname_H.mci.

home
fid = fopen(sprintf('%s%s_H.mci',directoryPath,myname),'r');
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
s(9).s = 'launchflag';
s(10).s = 'xs';
s(11).s = 'ys';
s(12).s = 'zs';
s(13).s = 'xfocus';
s(14).s = 'yfocus';
s(15).s = 'zfocus';
s(16).s = 'ux0';
s(17).s = 'uy0';
s(18).s = 'uz0';
s(19).s = 'radius';
s(20).s = 'waist';
s(21).s = 'Nt';

for i=1:21
    disp(sprintf('%d\t%10s = %0.4f',i,s(i).s,B(i)))
end

for j=1:B(21)
    i=i+1;
    disp(sprintf('---'))
    disp(sprintf('%d\tmua = %0.4f',i,B(i)))
    i=i+1;
    disp(sprintf('%d\tmus = %0.4f',i,B(i)))
    i=i+1;
    disp(sprintf('%d\tg   = %0.4f',i,B(i)))
end
