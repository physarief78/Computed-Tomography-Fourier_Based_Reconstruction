import numpy as np
import matplotlib.pyplot as plt

def load_matrix(filename, complex_data=False):
    """
    Load a matrix from a file.
    If complex_data is True, parse the entries as complex numbers.
    """
    with open(filename, 'r') as f:
        lines = f.readlines()

    if complex_data:
        # Parse complex numbers formatted as a+bi
        matrix = []
        for line in lines:
            row = []
            for entry in line.split():
                try:
                    # Replace "i" with "j" and parse as complex
                    entry = entry.replace("i", "j").replace(" ", "")
                    row.append(complex(entry))
                except ValueError:
                    print(f"Warning: Skipping malformed entry '{entry}'")
                    row.append(0)  # Add a placeholder for malformed data
            matrix.append(row)
    else:
        # Parse as floats
        matrix = [
            [float(entry) for entry in line.split()]
            for line in lines
        ]

    return np.array(matrix)

# Load the data
gaussian_signal = load_matrix("phantom_image.txt")
sinogram = load_matrix("sinogram.txt")
# fft_result = load_matrix("fourier_space.txt")
reconstructed_signal = load_matrix("reconstructed_image.txt")

# Visualize the Gaussian signal
plt.figure(figsize=(16, 12))

plt.subplot(1, 3, 1)
plt.title("Test Image")
plt.imshow(gaussian_signal, cmap='grey', interpolation='none')
# plt.colorbar(label="Intensity")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")

# # Visualize the FFT result (magnitude)
plt.subplot(1, 3, 2)
plt.title("Sinogram of Test Image")
plt.imshow(sinogram, cmap='grey', interpolation='none')
# plt.colorbar(label="Intensity")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")

#Visualize the reconstructed signal
plt.subplot(1, 3, 3)
plt.title("Reconstructed Signal")
plt.imshow(reconstructed_signal, cmap='grey', interpolation='none')
# plt.colorbar(label="Intensity")
plt.xlabel("X-axis")
plt.ylabel("Y-axis")

plt.tight_layout()
plt.show()
