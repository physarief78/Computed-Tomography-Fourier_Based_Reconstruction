# Computed-Tomography-Fourier_Based_Reconstruction
## 1. Introduction
Fourier based reconstruction is indeed fast and was one of the early methods used in CT image reconstruction. However, it has largely been superseded by more sophisticated algorithms such as FBP and iterative methods that handle real-world complexities better, including noise and artifacts. For certain applications where speed is critical and the data conditions are ideal, Fourier-based methods can still be very effective.

In this repository ihave done the Fourier-based reconstruction for magnetic particle imaging (MPI). This was done because MPI is a cutting-edge imaging technique with great potential in medical diagnostics and research. It allows for high-resolution, real-time imaging of magnetic nanoparticles, making it valuable for applications such as cancer detection, vascular imaging, for safety, and precision.

## 2. Mathematical Expression:
### 2.1. Signal Acquisition
The mathematical expression that be used in this method is
$$D'[x][y] = D[x][y] \cdot e^{i \phi(x,y)}$$
where x and y is the index that represent 2D matrix. Then, the defintion of $\phi$ is
$$\phi(x,y)=2 \pi f \frac{x + y}{N_r}$$
where $N_r$ is the total number of rows. If we assume there are some noises in the signal acquisition process, then we applied the noise function. Consider the noise is matching with Gaussian distribution, the mathematical expression for this noise is
$$G(x,y) = N(\mu, \sigma) + i N(\mu, \sigma)$$
N is the representation of a sample from Gaussian distribution with mean $\mu$ and the standard deviasion $\sigma$.
Then, the final mathematical expression of signal acquisition is
$$D''[x][y] = D'[x][y] + G(x,y)$$

### 2.2. Cooley-Tukey Radix-2 Decimation-in-Time (DIT) Fast Fourier Transform (FFT) Algorithm
For implement the Cooley-Tukey algortihm, we must define the \textit{discrete fourier transform} DFT algortihm first. The DFT has mathematical expression of
$$X[m] = \sum_{k=0}^{n-1} x[k] \cdot e^{-2\pi i \frac{mk}{n}}, \quad m = 0, 1, \dots, n-1$$
the inverse mathematical expression of this is
$$x[k] = \frac{1}{n} \sum_{m=0}^{n-1} X[m] \cdot e^{2\pi i \frac{mk}{n}}, \quad k = 0, 1, \dots, n-1$$

For doing faster computation, after defined the DFT expression, now we implement the Cooley-Tukey algorithm.
1. **Devide Step**
- Split the input \( x[k] \) into two smaller arrays:
     - $$ x_{\text{even}}[k] = x[2k] $$: Elements at even indices.
     - $$ x_{\text{odd}}[k] = x[2k+1] $$: Elements at odd indices.
   - If the size of $$ x[k] $$ is $$ n $$, the size of both $$ x_{\text{even}}[k] $$ and $$ x_{\text{odd}}[k] $$ is $$ n/2 $$.

2. **Recursive FFT**
 - Compute the FFT of the even and odd parts:
     
     $$X_{\text{even}}[m] = \text{FFT}(x_{\text{even}}[k])$$
     
     $$X_{\text{odd}}[m] = \text{FFT}(x_{\text{odd}}[k])$$
   
3. **Combine Step**:
   - Combine \( X_{\text{even}}[m] \) and \( X_{\text{odd}}[m] \) using the symmetry and periodicity properties of the DFT:
     \[
     X[m] = X_{\text{even}}[m] + W_m \cdot X_{\text{odd}}[m]
     \]
     \[
     X[m + n/2] = X_{\text{even}}[m] - W_m \cdot X_{\text{odd}}[m]
     \]
     - Here, \( W_m = e^{-2\pi i \frac{m}{n}} \) is the complex "twiddle factor".
   - For the **inverse FFT**, the twiddle factor becomes:
     \[
     W_m = e^{2\pi i \frac{m}{n}}
     \]
     - The results are also divided by 2 at each level of recursion for normalization.

4. **Base Case**:
   - If the size of the input \( n = 1 \), the FFT result is simply:
     \[
     X[0] = x[0]
     \]
