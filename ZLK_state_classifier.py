import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

from scipy.signal import butter, filtfilt

# Get current working directory
current_directory = os.getcwd()
print("Current working directory:", current_directory)
# Extract the last part of the path
last_part = os.path.basename(current_directory)

# Convert it to an integer
try:
    last_part_int = int(last_part)
    print("Extracted integer:", last_part_int)
except ValueError:
    print("The last part of the directory is not an integer.")
    
nin=(last_part_int-1)*3

print(nin)

# Set the sampling frequency and a lower cutoff frequency 
fs = 1/600  
fc = 0.000001   

# Scale fc to equivalent Nyquist units
nyquist_freq = 0.5 * fs   
cutoff = fc / nyquist_freq  
b, a = butter(4, cutoff, btype='lowpass')

# Here it starts the main loop

with open('results.txt', 'w') as f:
    for i in range(1,4):
        j = nin + i
        filename1 = 'res_arg_{:02d}'.format(i)
        filename2 = 'res_fil_{:07d}.png'.format(i)
        filename3 = 'res_osc_{:07d}.png'.format(i)
        print(filename1,filename2,filename3)
    
        # Load your signal data into a NumPy array
        df = pd.read_csv(filename1,
                         skiprows=0,
                         header=None,
                         sep = '\s+',
                         index_col=None,
                         low_memory=False,
                         names=['time','cos','sin','angle'],
                         dtype={'time':np.float64,
                                'cos':np.float64,
                                'sin':np.float64,
                                'angle':np.float64
                                })
        data= df['angle'].values
        # Apply a Butterworth low-pass filter to extract the long-period frequencies  
        filtered_data = filtfilt(b, a, data)
        t = np.arange(0, len(data)) / fs
    
        # Plot the filtered resonant angle
        plt.plot(t, filtered_data, label='Filtered Signal')
        plt.xlabel('Time [yr]')
        plt.ylabel('Filtered resonant angle [degrees]')
        plt.ylim([0,360])
        plt.savefig(filename2)
        plt.close()
        # Plot the osculating resonant angle
        plt.plot(t, data, label='Original Signal')
        plt.xlabel('Time [yr]')
        plt.ylabel('Osculating Resonant angle [degrees]')
        plt.savefig(filename3)
        plt.close()
        delta = max(data) - min(data)
        mean_a = np.mean(data)-180.0
        if mean_a < 0.0:
            mean_a = mean_a +360.0
        
        # Calculate the power spectrum
        fft_data = np.fft.fft(filtered_data)
        power_spectrum = np.abs(fft_data) ** 2
        n = data.size
        freq = np.fft.fftfreq(n) * fs
    
        # Find the indices of the largest amplitudes in descending order
        sorted_indices = np.argsort(power_spectrum[:n//2])[::-1]
        dominant_freq = freq[sorted_indices[0]]
        if dominant_freq != 0 :
            dominant_period = 1 / dominant_freq
            f.write('{:02d} {:.2f} {:.2f} {:.2f}\n'.format(j, mean_a, delta, dominant_period/1e6))
        else :
        
            second_max_amp_idx = sorted_indices[1]  # Index of the second largest amplitudefreq = freq[freq != 0]
    
            # Calculate the second dominant frequency and period
            second_dominant_freq = freq[second_max_amp_idx]
            second_dominant_period = 1 / second_dominant_freq
            f.write('{:02d} {:.2f} {:.2f} {:.2f}\n'.format(j, mean_a, delta, second_dominant_period/1e6))

