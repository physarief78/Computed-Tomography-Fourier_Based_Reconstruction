import numpy as np
import matplotlib.pyplot as plt

# Load data from the text files
simulated_signal = np.loadtxt('phantom_image.txt')
withoutfilter_signal = np.loadtxt('reconstructed_withoutFilter_image.txt')
reconstructed_image = np.loadtxt('reconstructed_image.txt')

# Plot simulated signal
plt.figure(figsize=(12, 4))
plt.subplot(1,3,1)
plt.title('Original Phantom')
plt.xlabel('Rows')
plt.ylabel('Columns')
plt.legend()
plt.imshow(simulated_signal, cmap='gray')

# Plot deconvolved signal
plt.subplot(1, 3, 2)
plt.title('Reconstructed Image Without Filter')
plt.xlabel('Rows')
plt.ylabel('Colmuns')
plt.legend()
plt.imshow(withoutfilter_signal, cmap='gray')

# Display reconstructed image
plt.subplot(1, 3, 3)
plt.imshow(reconstructed_image, cmap='gray')
plt.title('Reconstructed Image')
plt.show()
