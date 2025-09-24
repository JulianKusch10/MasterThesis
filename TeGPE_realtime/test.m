N = 50;
h = waitbar(0,'Please wait...');
for k = 1:N
    pause(0.05);
    waitbar(k/N,h,sprintf('Step %d of %d',k,N));
end
close(h)