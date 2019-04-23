
%% Test to determine optimal block size
% A = rand(4*3000,4*3000,'single');
% 
% tic
% for j=1:10
% 	B = A.';
% end
% time_matlab = toc;
% 
% n = 13;
% speedups = zeros(1,n);
% for i=1:n
% 	4*2^(i-1)
% 
% 	tic
% 	for j=1:10
% 		C = homemadepermute(A,4*2^(i-1),true);
% 	end
% 	time_C = toc;
% 
% 	speedups(i) = 100*(time_matlab/time_C-1);
% end
% 
% semilogx(4*2.^((1:n)-1),speedups)

%% Test to determine speedup as function of matrix size
n = 45;
N = round(logspace(log10(10),log10(3000),n));

speedups = zeros(1,n);
for i=1:n
	4*N(i)
	A = rand(4*N(i),4*N(i),'single');

	tic
	for j=1:10
		B = homemadepermute(A,256,true);
	end
	time_SSE = toc;

	tic
	for j=1:10
		C = homemadepermute(A,256,false);
	end
	time_noSSE = toc;

	speedups(i) = time_SSE/time_noSSE;
end

semilogx(4*N,speedups)

% 	fprintf('Speedup is %.1f%%\n',100*(time_matlab/time_C-1));

% A = rand(8,8,'single')
% B = homemadepermute(A,256,false)