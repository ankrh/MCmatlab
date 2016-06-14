%Program to input graph files in a folder, made with the ocean optics
%spectrometer.  
clc;
clear all

ext = '.txt'; % Only find files with this extention
AllFiles = dirrec(char(fileparts(pwd)),ext); % dirrec finds all files in the folder and it's sub folder with the given extention
name = 'VL2_9,2J';
name2 = 'VL+_9,3J';

k=0;
for i = 1:length(AllFiles)
    if isempty(strfind(char(AllFiles(i)),[name]))==0
        k = k+1;
        files(k)=AllFiles(i);
    elseif isempty(strfind(char(AllFiles(i)),[name2]))==0
        k = k+1;
        files(k)=AllFiles(i);
    end
end


for n = 1:k
    fid = fopen(char(files(n)),'r');
    A = textscan(fid, '%8f %8f', 'HeaderLines', 0,'Delimiter','\t');
    fclose(fid);
    clear fid

    if n == 1 
        nm1 = A{1};
        c1 = A{2};
    elseif n == 2
        nm2 = A{1};
        c2 = A{2};
    end

end

M = max(c1)-c1(1); %max of VL+ minus the background

a = 0;
b = 0;
c = 0;
d = 0;
e = 0;
f = 0;

cc1 = zeros(length(nm1),1);
cc2 = zeros(length(nm1),1);
cc3 = zeros(length(nm1),1);

for n = 1:length(nm1)
    if nm1(n) <= 500 
        a = a+1;
        cc1(n) = 0;
    elseif nm1(n) > 500 && nm1(n) < 600
        b = b+1;
        cc1(n) = (n-a)*(M/9);
    elseif nm1(n) >=600 && nm1(n) <= 800
        cc1(n) = M - (n-a-b)*(M/21);
    elseif nm1(n) > 800
        cc1(n) = 0;
    end
    
    if nm1(n) <= 490 
        c = c+1;
        cc2(n) = 0;
    elseif nm1(n) > 490 && nm1(n) < 590
        d = d+1;
        cc2(n) = (n-c)*(M/9);
    elseif nm1(n) == 590
        d = d+1;
        cc2(n) = cc2(n-1);
    elseif nm1(n) >=600 && nm1(n) <= 800
        cc2(n) = M - (n-c-d)*(M/21);
    elseif nm1(n) > 800
        cc2(n) = 0;
    end
    
    if nm1(n) <= 500 
        e = e+1;
        cc3(n) = 0;
    elseif nm1(n) > 500 && nm1(n) < 600
        f = f+1;
        cc3(n) = (n-e)*(M/9);
    elseif nm1(n) == 600
        f = f+1;
        cc3(n) = cc3(n-1);
    elseif nm1(n) >=610 && nm1(n) <= 810
        cc3(n) = M - (n-e-f)*(M/21);
    elseif nm1(n) > 810
        cc3(n) = 0;
    end
end
    
% bg1 = sum(c1(1:5)/5);
% bg2 = sum(c2(1:5)/5);
% 
% cc1 = c1 - bg1;
% cc2 = c2 - bg2;

%%
Total1 = sum(cc1); %mid
Total2 = sum(cc2); %left
Total3 = sum(cc3); %right

norm2 = Total1/Total2;
norm3 = Total1/Total3;

counts1 = cc1;
counts2 = cc2.*norm2;
counts3 = cc3.*norm3;

% save('test2','nm','Power')
% save(['VL+_norm_' p 'J'],'nm2','counts2')
%%
fid1 = fopen('Tri_mid.txt','wt');
A = [nm1'; counts1'];
fprintf(fid1,'%7.2f %12.4f\n', A); 
fclose(fid1);

fid2 = fopen('Tri_left.txt','wt');
B = [nm1'; counts2'];
fprintf(fid2,'%7.2f %11.4f\n', B); 
fclose(fid2);

fid3 = fopen('Tri_right.txt','wt');
C = [nm1'; counts3'];
fprintf(fid3,'%7.2f %11.4f\n', C); 
fclose(fid3);
    
 



