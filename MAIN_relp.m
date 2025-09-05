% This is a simulation of the Regular Pulse Excited Long Term Prediction
% algorithm used in GSM 6.10.

%MAIN BODY

close all;
clear all;
clc;

inpfilenm = 't_8k_2s';
[x, fs] =wavread(inpfilenm); 

%or record your own voice files using,
% t = 2; 		% = how many seconds
% fs = 8000;    % = sampling frequency
% x = wavrecord(t.*fs, fs, 1);

%LENGTH DISPLAY (IN SEC) OF INPUT WAVEFILE,
t=length(x)./fs;
sprintf('Processing the wavefile "%s"', inpfilenm)
sprintf('The wavefile is  %3.2f  seconds long', t)


%COMPRESSION STARTS HERE,
[aCoeff, b_LTopt, Topt, e_prime] = f_ENCODER_relp(x, fs);
        % e_prime is instead of position, peak_magitude_index and sample_amplitude_index. (temporarily)

[synth_speech] = f_DECODER_relp(aCoeff, b_LTopt, Topt, e_prime);




%RESULTS,
beep;
disp('Press a key to play the original sound!');
pause;
soundsc(x, fs);

disp('Press a key to play the RPE-LTP (RELP) compressed sound!');
pause;
soundsc(synth_speech, fs);

figure;
subplot(211),plot(x); title(['Original signal = "', inpfilenm, '"']);
subplot(212), plot(synth_speech); title('RELP compressed output');
