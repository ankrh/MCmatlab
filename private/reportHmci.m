function H_mci = reportHmci(directoryPath,myname)
% function reportHmci(directoryPath,myname)
%   Lists the values of the input file  myname_H.mci.
fprintf('------ mcxyz %s -------\n',myname)

fid = fopen(sprintf('%s%s_H.mci',directoryPath,myname),'r');
B = fscanf(fid,'%f');
fclose(fid);

H_mci.time_min      = B(1);
H_mci.nx            = B(2);
H_mci.ny            = B(3);
H_mci.nz            = B(4);
H_mci.dx            = B(5);
H_mci.dy            = B(6);
H_mci.dz            = B(7);
H_mci.beamtypeflag  = B(8);
H_mci.boundaryflag  = B(9);
H_mci.xfocus        = B(10);
H_mci.yfocus        = B(11);
H_mci.zfocus        = B(12);
H_mci.ux0           = B(13);
H_mci.uy0           = B(14);
H_mci.uz0           = B(15);
H_mci.waist         = B(16);
H_mci.divergence    = B(17);
H_mci.Nt            = B(18);

names = fieldnames(H_mci);
for i=1:18
    fprintf('%d\t%10s = %0.8f\n',i,names{i},B(i))
end

for j=1:B(18)
    i=i+1;
    fprintf('---\n')
    fprintf('%d\tmua = %0.8f\n',i,B(i))
    Hmci.mua(j) = B(i);
    i=i+1;
    fprintf('%d\tmus = %0.8f\n',i,B(i))
    Hmci.mus(j) = B(i);
    i=i+1;
    fprintf('%d\tg   = %0.8f\n',i,B(i))
    Hmci.g(j) = B(i);
end
