function [c_eff, output_amp, lengths_vec, wavenumber_vec] = fft_lengths(amplitude_vec, freq_vec, mu_eff, epsilon_eff)

% the function recieves a frequenct domain map (frq_vac amp_vec) and
% computes it's fft and then using the provided mu and epsilon chnages it
% to the length domain.

start_freq = freq_vec(1); end_freq = freq_vec(end); L = length(freq_vec); 
Fs = freq_vec(2) - freq_vec(1);
output_amp = fft(amplitude_vec);

char_time_vec = linspace(0, 1/Fs, L);
c_eff = (mu_eff * epsilon_eff)^(-0.5);

lengths_vec = char_time_vec * c_eff;
wavenumber_vec = freq_vec / c_eff;

end

