function H_mci = reportHmci(directoryPath,myname)
% function reportHmci(directoryPath,myname)
%   Lists the values of the input file  myname_H.mci.

fid = fopen(sprintf('%s%s_H.mci',directoryPath,myname),'r');
B = fscanf(fid,'%f');
fclose(fid);

H_mci(1).name = 'time_min';
H_mci(2).name = 'Nx';
H_mci(3).name = 'Ny';
H_mci(4).name = 'Nz';
H_mci(5).name = 'dx';
H_mci(6).name = 'dy';
H_mci(7).name = 'dz';
H_mci(8).name = 'mcflag';
H_mci(9).name = 'launchflag';
H_mci(10).name = 'xs';
H_mci(11).name = 'ys';
H_mci(12).name = 'zs';
H_mci(13).name = 'xfocus';
H_mci(14).name = 'yfocus';
H_mci(15).name = 'zfocus';
H_mci(16).name = 'ux0';
H_mci(17).name = 'uy0';
H_mci(18).name = 'uz0';
H_mci(19).name = 'radius';
H_mci(20).name = 'waist';
H_mci(21).name = 'Nt';

for i=1:21
    H_mci(i).value = B(i);
    fprintf('%d\t%10s = %0.4f\n',i,H_mci(i).name,B(i))
end

for j=1:B(21)
    i=i+1;
    fprintf('---\n')
    fprintf('%d\tmua = %0.4f\n',i,B(i))
    i=i+1;
    fprintf('%d\tmus = %0.4f\n',i,B(i))
    i=i+1;
    fprintf('%d\tg   = %0.4f\n',i,B(i))
end
