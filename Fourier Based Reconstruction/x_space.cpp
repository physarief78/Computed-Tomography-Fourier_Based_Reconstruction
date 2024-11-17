#include <iostream>
#include <fstream>
#include <vector>
#include <complex>
#include <cmath>
#include <chrono>
#include <random>

using namespace std;

const double PI = 3.14159265358979323846;

// Function to load Shepp-Logan phantom data
vector<vector<double>> loadPhantom(const string &filename, int &rows, int &cols) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }

    rows = 512; // Assuming size based on the problem description
    cols = 512;

    vector<vector<double>> phantom(rows, vector<double>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (!(file >> phantom[i][j])) {
                cerr << "Error: Missing or invalid data at row " << i << ", column " << j << endl;
                exit(1);
            }
        }
    }

    file.close();
    cout << "Successfully loaded phantom image with dimensions: " << rows << "x" << cols << endl;
    return phantom;
}

// Function to add Gaussian noise
void addGaussianNoise(vector<vector<complex<double>>> &data, double mean, double stddev) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<double> dist(mean, stddev);

    int rows = data.size();
    int cols = data[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double noiseReal = dist(gen);
            double noiseImag = dist(gen);
            data[i][j] += complex<double>(noiseReal, noiseImag);
        }
    }

    cout << "Added Gaussian noise with mean: " << mean << " and stddev: " << stddev << endl;
}

// Apply magnetic field effect with added noise
void applyMagneticFieldWithNoise(vector<vector<complex<double>>> &data, double frequency, double noiseMean, double noiseStddev) {
    int rows = data.size();
    int cols = data[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double phase = 2 * PI * frequency * (i + j) / rows;
            complex<double> phaseShift(cos(phase), sin(phase));
            data[i][j] *= phaseShift;
        }
    }

    // Add noise after applying the magnetic field
    addGaussianNoise(data, noiseMean, noiseStddev);
}

// Undo magnetic field effect (reverse sinusoidal phase shift)
void undoMagneticField(vector<vector<complex<double>>> &data, double frequency) {
    int rows = data.size();
    int cols = data[0].size();

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double phase = 2 * PI * frequency * (i + j) / rows;
            complex<double> phaseShift(cos(-phase), sin(-phase));
            data[i][j] *= phaseShift;
        }
    }
}

// Perform 1D FFT on a vector
void fft1D(vector<complex<double>> &data, bool inverse = false) {
    int n = data.size();
    if ((n & (n - 1)) != 0) {
        cerr << "Error: FFT input size must be a power of 2. Got size " << n << endl;
        exit(1);
    }

    if (n <= 1) return;

    vector<complex<double>> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i) {
        even[i] = data[i * 2];
        odd[i] = data[i * 2 + 1];
    }

    fft1D(even, inverse);
    fft1D(odd, inverse);

    double angle = 2 * PI / n * (inverse ? 1 : -1);
    complex<double> w(1), wn(cos(angle), sin(angle));
    for (int i = 0; i < n / 2; ++i) {
        data[i] = even[i] + w * odd[i];
        data[i + n / 2] = even[i] - w * odd[i];
        if (inverse) {
            data[i] /= 2;
            data[i + n / 2] /= 2;
        }
        w *= wn;
    }
}

// FFT 2D
template <typename T>
vector<vector<complex<double>>> fft2D(const vector<vector<T>> &image, bool inverse = false) {
    int rows = image.size();
    int cols = image[0].size();

    vector<vector<complex<double>>> transformed(rows, vector<complex<double>>(cols));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            transformed[i][j] = static_cast<complex<double>>(image[i][j]);
        }
    }

    // FFT rows
    for (int i = 0; i < rows; ++i)
        fft1D(transformed[i], inverse);

    // FFT columns
    for (int j = 0; j < cols; ++j) {
        vector<complex<double>> column(rows);
        for (int i = 0; i < rows; ++i)
            column[i] = transformed[i][j];

        fft1D(column, inverse);

        for (int i = 0; i < rows; ++i)
            transformed[i][j] = column[i];
    }

    return transformed;
}

// Convert complex matrix to real
vector<vector<double>> complexToReal(const vector<vector<complex<double>>> &data) {
    int rows = data.size();
    int cols = data[0].size();

    vector<vector<double>> realData(rows, vector<double>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            realData[i][j] = data[i][j].real();

    return realData;
}

vector<vector<double>> wienerFilter(const vector<vector<double>> &image, int windowSize) {
    int rows = image.size();
    int cols = image[0].size();
    int halfWindow = windowSize / 2;

    vector<vector<double>> filtered(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            double localMean = 0.0;
            double localVariance = 0.0;
            double pixelCount = 0.0;

            // Compute local mean and variance
            for (int wi = -halfWindow; wi <= halfWindow; ++wi) {
                for (int wj = -halfWindow; wj <= halfWindow; ++wj) {
                    int ni = i + wi;
                    int nj = j + wj;
                    if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
                        localMean += image[ni][nj];
                        ++pixelCount;
                    }
                }
            }
            localMean /= pixelCount;

            for (int wi = -halfWindow; wi <= halfWindow; ++wi) {
                for (int wj = -halfWindow; wj <= halfWindow; ++wj) {
                    int ni = i + wi;
                    int nj = j + wj;
                    if (ni >= 0 && ni < rows && nj >= 0 && nj < cols) {
                        localVariance += pow(image[ni][nj] - localMean, 2);
                    }
                }
            }
            localVariance /= pixelCount;

            double noiseVariance = 1e1; // Adjust based on expected noise level
            double pixelValue = image[i][j];
            double restoredPixel = localMean + max(0.0, 1.0 - noiseVariance / localVariance) * (pixelValue - localMean);

            filtered[i][j] = restoredPixel;
        }
    }

    cout << "Applied Wiener filter with window size: " << windowSize << endl;
    return filtered;
}

vector<vector<double>> totalVariationDenoising(const vector<vector<double>> &image, double lambda, int iterations) {
    int rows = image.size();
    int cols = image[0].size();
    vector<vector<double>> denoised = image; // Initialize with the noisy image

    // Gradient descent for TV minimization
    for (int iter = 0; iter < iterations; ++iter) {
        vector<vector<double>> grad(rows, vector<double>(cols, 0.0));

        // Compute gradient of TV term
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                double gradX = (i + 1 < rows) ? denoised[i + 1][j] - denoised[i][j] : 0.0;
                double gradY = (j + 1 < cols) ? denoised[i][j + 1] - denoised[i][j] : 0.0;

                double norm = sqrt(gradX * gradX + gradY * gradY + 1e-8); // Avoid division by zero
                grad[i][j] += (gradX / norm) + (gradY / norm);
            }
        }

        // Update the image using the gradient
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                denoised[i][j] -= lambda * grad[i][j]; // Gradient descent update
            }
        }
    }

    cout << "Applied Total Variation Denoising with lambda: " << lambda << " and iterations: " << iterations << endl;
    return denoised;
}

// Save image to file
void saveImageToText(const string &filename, const vector<vector<double>> &image) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open file " << filename << endl;
        exit(1);
    }

    for (const auto &row : image) {
        for (double pixel : row) {
            file << pixel << " ";
        }
        file << endl;
    }

    file.close();
    cout << "Image saved to " << filename << endl;
}

int main() {
    int rows, cols;
    string inputFile = "phantom_image.txt";

    auto start = chrono::high_resolution_clock::now();

    // Step 1: Load phantom data
    vector<vector<double>> phantom = loadPhantom(inputFile, rows, cols);

    // Convert to complex for further processing
    vector<vector<complex<double>>> phantomComplex(rows, vector<complex<double>>(cols));
    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            phantomComplex[i][j] = complex<double>(phantom[i][j]);

    // Step 2: Apply magnetic field with noise
    double frequency = 1000.0; // Example frequency for simulation
    double noiseMean = 0.01;    // Noise mean
    double noiseStddev = 0.25; // Noise standard deviation
    applyMagneticFieldWithNoise(phantomComplex, frequency, noiseMean, noiseStddev);

    // Step 3: Undo magnetic field during reconstruction
    undoMagneticField(phantomComplex, frequency);

    // Step 4: FFT of the acquired phantom data from the magnetic field acquitition
    vector<vector<complex<double>>> fftResult = fft2D(phantomComplex);

    // Step 5: Inverse FFT to reconstruct the signal
    vector<vector<complex<double>>> reconstructedComplex = fft2D(fftResult, true);
    vector<vector<double>> reconstructed = complexToReal(reconstructedComplex);

    string outputFile1 = "reconstructed_withoutFilter_image.txt";
    saveImageToText(outputFile1, reconstructed);

    // Step 6: Apply Wiener filter for noise reduction
    int windowSize = 7.0; // Size of the local window
    reconstructed = wienerFilter(reconstructed, windowSize);

    double lambda = 1e-15;    // Regularization parameter
    int iterations = 100;   // Number of gradient descent iterations
    reconstructed = totalVariationDenoising(reconstructed, lambda, iterations);

    // Step 7: Save reconstructed image
    string outputFile2 = "reconstructed_image.txt";
    saveImageToText(outputFile2, reconstructed);

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "Total execution time: " << elapsed.count() << " seconds" << endl;

    return 0;
}


