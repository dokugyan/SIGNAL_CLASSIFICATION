function [signalRecon] = stftRecon(S, Fs, win, nWindows, overlap)

[signalRecon,T] = istft(S,Fs,'Window',win,'OverlapLength',overlap);
        % Remove zeros.
signalRecon(1:nWindows) = [];
signalRecon(end-nWindows+1:end) = [];
T = T(1:end-2*nWindows);


