from scipy import signal as sig
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt


# Generate Sin Wave
L, f, f_s = 1000, 50, 2e3
t = np.arange(0, L)
y = np.sin(2*math.pi*f*t/f_s)

plt.plot(t, y)
plt.show()
# print y

print('Sithan Made a Change at 10pm')
print('Antonello Made a Change on master at 10pm')

# Ignore this: Just creates a signal

L, f, f_s = 1000, 50, 2e3
t = np.arange(0, L)
x = np.sin(2*math.pi*f*t/f_s)
plt.plot(x)
plt.show()

# Fourier Transfrom is called FFT : It stands for Fast Fourier Transfrom
X = np.fft.fft(x)

# You take the absolute value of the Fourier Transform since the FFT gives you a complex number
plt.stem(np.abs(X))

# The second half of the FFT signal is a mirror image of the first half and can be igonored
plt.xlim([0, 500])

plt.show()

# There is only one spike because your original signal has only one sine wave.
# If you increase the number of sine wave you will have more spikes. See this example with 2 sine waves:

x_new = x + 0.5*np.sin(2*math.pi*(f+100)*t/f_s)
X_new = np.fft.fft(x_new)
plt.stem(np.abs(X_new))
plt.show()
