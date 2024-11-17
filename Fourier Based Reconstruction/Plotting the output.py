# import numpy as np
# import matplotlib.pyplot as plt

# # Load the phantom image
# phantom = np.loadtxt("simulated_signal.txt")

# #
# signal = np.loadtxt("deconvolved_signal.txt")

# # Load the reconstructed image
# reconstructed_image = np.loadtxt("reconstructed_image.txt")

# # Simulate the detected signal by convolving the phantom with the PSF
# # Here we approximate by reading the simulated signal directly, if available.
# # Uncomment and modify the following if you have generated the intermediate signal file in C++.
# # signal = np.loadtxt("signal_image.txt")

# # Display the images
# plt.figure(figsize=(15, 5))

# # Original Phantom Image
# plt.subplot(1, 3, 1)
# plt.title("Original Phantom Image")
# plt.imshow(phantom, cmap='gray')
# plt.colorbar()

# # Uncomment below if you have simulated `signal_image.txt` from C++:
# # Simulated Signal
# plt.subplot(1, 3, 2)
# plt.title("Simulated Signal (Convolved)")
# plt.imshow(signal, cmap='gray')
# plt.colorbar()

# # Reconstructed Image
# plt.subplot(1, 3, 3)
# plt.title("Reconstructed Image")
# plt.imshow(reconstructed_image, cmap='gray')
# plt.colorbar()

# plt.tight_layout()
# plt.show()

# import numpy as np
# import matplotlib.pyplot as plt

# # Load the normalized 2D signal and reconstructed image from text files
# signal_2d = np.loadtxt("signal.txt")
# reconstructed_image_2d = np.loadtxt("reconstructed_image.txt")

# # Plot the normalized 2D signal
# plt.figure(figsize=(12, 6))
# plt.subplot(1, 2, 1)
# plt.imshow(signal_2d, cmap="viridis", origin="lower")
# plt.colorbar()
# plt.title("Normalized 2D Magnetic Response Signal")
# plt.xlabel("X Position")
# plt.ylabel("Y Position")

# # Plot the normalized 2D reconstructed image
# plt.subplot(1, 2, 2)
# plt.imshow(reconstructed_image_2d, cmap="viridis", origin="lower")
# plt.colorbar()
# plt.title("Normalized 2D Reconstructed Image")
# plt.xlabel("X Position")
# plt.ylabel("Y Position")

# # Display the plots
# plt.tight_layout()
# plt.show()

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
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.legend()
plt.imshow(simulated_signal, cmap='gray')

# Plot deconvolved signal
plt.subplot(1, 3, 2)
plt.title('Reconstructed Image Without Filter')
plt.xlabel('Sample Index')
plt.ylabel('Amplitude')
plt.legend()
plt.imshow(withoutfilter_signal, cmap='gray')

# Display reconstructed image
plt.subplot(1, 3, 3)
plt.imshow(reconstructed_image, cmap='gray')
plt.title('Reconstructed Image')
plt.show()