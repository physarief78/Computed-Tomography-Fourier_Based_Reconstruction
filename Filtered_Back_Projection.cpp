#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <omp.h>
#include <complex>
#include <chrono>
#include <algorithm>

constexpr double PI = 3.141592653589793;
using Complex = std::complex<double>;
using Clock = std::chrono::high_resolution_clock;
using VecComplex = std::vector<Complex>;

// Function to load data from a .txt file
std::vector<std::vector<double>> loadData(const std::string& filename) {
    std::ifstream inputFile(filename);
    if (!inputFile) {
        std::cerr << "Error opening file.\n";
        exit(1);
    }

    std::vector<std::vector<double>> data;
    std::string line;

    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        std::vector<double> row;
        double value;

        while (ss >> value) {
            row.push_back(value);
        }

        if (!row.empty()) {
            data.push_back(row);
        }
    }

    inputFile.close();
    return data;
}

// FFT Implementation
void bitReversePermutation(VecComplex &data) {
    int n = data.size();
    int log2n = log2(n);

    for (int i = 0; i < n; ++i) {
        int reversed = 0;
        for (int j = 0; j < log2n; ++j) {
            if (i & (1 << j)) {
                reversed |= (1 << (log2n - 1 - j));
            }
        }
        if (i < reversed) {
            std::swap(data[i], data[reversed]);
        }
    }
}

void fft(VecComplex &data, bool inverse = false) {
    int n = data.size();
    bitReversePermutation(data);

    for (int s = 1; s <= log2(n); ++s) {
        int m = 1 << s;
        Complex wm = std::polar(1.0, (inverse ? 2 : -2) * PI / m);
        for (int k = 0; k < n; k += m) {
            Complex w = 1.0;
            for (int j = 0; j < m / 2; ++j) {
                Complex t = w * data[k + j + m / 2];
                Complex u = data[k + j];
                data[k + j] = u + t;
                data[k + j + m / 2] = u - t;
                w *= wm;
            }
        }
    }

    if (inverse) {
        for (auto &x : data) {
            x /= n;
        }
    }
}

// Precompute sine and cosine lookup tables
void precomputeTrigTables(int numAngles, std::vector<double>& cosTable, std::vector<double>& sinTable) {
    cosTable.resize(numAngles);
    sinTable.resize(numAngles);

    for (int thetaIdx = 0; thetaIdx < numAngles; ++thetaIdx) {
        double theta = thetaIdx * PI / numAngles; // Convert angle to radians
        cosTable[thetaIdx] = std::cos(theta);
        sinTable[thetaIdx] = std::sin(theta);
    }
}

// Radon transform with precomputed trigonometric lookup
std::vector<std::vector<double>> radon_transform(const std::vector<std::vector<double>>& image, 
                                                 int numAngles, 
                                                 const std::vector<double>& cosTable, 
                                                 const std::vector<double>& sinTable) {
    int size = image.size();
    int center = size / 2;
    std::vector<std::vector<double>> sinogram(numAngles, std::vector<double>(size, 0.0));

    #pragma omp parallel for
    for (int i = 0; i < numAngles; ++i) {
        double cosTheta = cosTable[i];
        double sinTheta = sinTable[i];
        for (int t = 0; t < size; ++t) {
            double sum_val = 0.0;
            for (int s = 0; s < size; ++s) {
                int x = static_cast<int>(center + (t - center) * cosTheta - (s - center) * sinTheta);
                int y = static_cast<int>(center + (t - center) * sinTheta + (s - center) * cosTheta);
                if (x >= 0 && x < size && y >= 0 && y < size) {
                    sum_val += image[y][x];
                }
            }
            sinogram[i][t] = sum_val;
        }
    }
    return sinogram;
}

// Reconstruction with interpolation
std::vector<std::vector<double>> filtered_back_projection(const std::vector<std::vector<double>> &sinogram, int size, int numAngles) {
    int center = size / 2;
    std::vector<std::vector<double>> reconstructed(size, std::vector<double>(size, 0.0));

    // Apply ramp filter
    std::vector<std::vector<double>> filteredSinogram = sinogram;
    #pragma omp parallel for
    for (int i = 0; i < numAngles; ++i) {
        int n = filteredSinogram[i].size();
        VecComplex projection(n);
        for (int j = 0; j < n; ++j) projection[j] = Complex(filteredSinogram[i][j], 0.0);
        fft(projection);

        for (int j = 0; j < n / 2; ++j) {
            double ramp = j * 2.0 / n;
            projection[j] *= ramp;
            projection[n - j - 1] *= ramp;
        }
        fft(projection, true);

        for (int j = 0; j < n; ++j) filteredSinogram[i][j] = projection[j].real();
    }

    // Back projection with interpolation
    #pragma omp parallel for
    for (int i = 0; i < numAngles; ++i) {
        double theta = i * PI / numAngles;
        double cosTheta = std::cos(theta);
        double sinTheta = std::sin(theta);

        for (int y = 0; y < size; ++y) {
            for (int x = 0; x < size; ++x) {
                double t = (x - center) * cosTheta + (y - center) * sinTheta;
                int tIdx = std::floor(t + center);
                double w = t + center - tIdx;

                if (tIdx >= 0 && tIdx + 1 < size) {
                    reconstructed[y][x] += (1 - w) * filteredSinogram[i][tIdx] + w * filteredSinogram[i][tIdx + 1];
                }
            }
        }
    }

    // Normalize
    double normFactor = PI / numAngles;
    for (auto &row : reconstructed) {
        for (auto &pixel : row) {
            pixel *= normFactor;
        }
    }

    return reconstructed;
}

// Function to save the sinogram or image data to a file
void saveData(const std::string& filename, const std::vector<std::vector<double>>& data) {
    std::ofstream outputFile(filename);
    if (!outputFile) {
        std::cerr << "Error creating file.\n";
        exit(1);
    }

    for (const auto& row : data) {
        for (size_t j = 0; j < row.size(); ++j) {
            outputFile << row[j];
            if (j < row.size() - 1) {
                outputFile << " ";
            }
        }
        outputFile << "\n";
    }

    outputFile.close();
}

int main() {
    // Load the input image
    std::string inputFilename = "phantom_image.txt";
    std::vector<std::vector<double>> image = loadData(inputFilename);

    // Precompute trigonometric tables
    int numAngles = 1080;
    std::vector<double> cosTable, sinTable;
    precomputeTrigTables(numAngles, cosTable, sinTable);

    // Compute the Radon transform
    auto start = Clock::now();
    std::vector<std::vector<double>> sinogram = radon_transform(image, numAngles, cosTable, sinTable);
    auto end = Clock::now();

    std::cout << "Radon transform completed in "
              << std::chrono::duration<double>(end - start).count()
              << " s.\n";

    // Save the sinogram
    std::string sinogramFilename = "sinogram.txt";
    saveData(sinogramFilename, sinogram);

    std::cout << "Sinogram saved to " << sinogramFilename << "\n";

    // Reconstruct the image using Filtered Back Projection
    start = Clock::now();
    std::vector<std::vector<double>> reconstructed = filtered_back_projection(sinogram, image.size(), numAngles);
    end = Clock::now();

    std::cout << "Reconstruction completed in "
              << std::chrono::duration<double>(end - start).count()
              << " s.\n";

    // Save the reconstructed image
    std::string reconstructedFilename = "reconstructed_image.txt";
    saveData(reconstructedFilename, reconstructed);

    std::cout << "Reconstructed image saved to " << reconstructedFilename << "\n";

    return 0;
}
