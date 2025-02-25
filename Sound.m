function[] = Sound(T)
Fs = 44100;  % Sampling frequency (Hz)
t = 0:1/Fs:T;  % Time vector (1 second long)
y = sin(2 * pi * 440 * t);  % 440 Hz tone (A4 note)
% Play sound
sound(y, Fs);