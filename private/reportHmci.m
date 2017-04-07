function H_mci = reportHmci(directoryPath,myname)
% function reportHmci(directoryPath,myname)
%   Lists the values of the input file  myname_H.mci.
fprintf('------ mcxyz %s -------\n',myname)

fid = fopen(sprintf('%s%s_H.mci',directoryPath,myname),'r');
B = fscanf(fid,'%f');
fclose(fid);

H_mci.time_min      = B(1);
H_mci.Nx            = B(2);
H_mci.Ny            = B(3);
H_mci.Nz            = B(4);
H_mci.dx            = B(5);
H_mci.dy            = B(6);
H_mci.dz            = B(7);
H_mci.mcflag        = B(8);
H_mci.launchflag    = B(9);
H_mci.boundaryflag  = B(10);
H_mci.xs            = B(11);
H_mci.ys            = B(12);
H_mci.zs            = B(13);
H_mci.xfocus        = B(14);
H_mci.yfocus        = B(15);
H_mci.zfocus        = B(16);
H_mci.ux0           = B(17);
H_mci.uy0           = B(18);
H_mci.uz0           = B(19);
H_mci.radius        = B(20);
H_mci.waist         = B(21);
H_mci.Nt            = B(22);

names = fieldnames(H_mci);
for i=1:22
    fprintf('%d\t%10s = %0.4f\n',i,names{i},B(i))
end

for j=1:B(22)
    i=i+1;
    fprintf('---\n')
    fprintf('%d\tmua = %0.4f\n',i,B(i))
    Hmci.mua(j) = B(i);
    i=i+1;
    fprintf('%d\tmus = %0.4f\n',i,B(i))
    Hmci.mus(j) = B(i);
    i=i+1;
    fprintf('%d\tg   = %0.4f\n',i,B(i))
    Hmci.g(j) = B(i);
end
