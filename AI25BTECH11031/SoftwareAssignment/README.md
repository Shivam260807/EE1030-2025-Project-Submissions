Image Compression Using Truncated SVD
Project Overview
This project implements image compression using the Truncated Singular Value Decomposition (SVD).
A grayscale image is represented as a matrix, where each element represents a pixel intensity (0 = black, 255 = white).
Using SVD, we can decompose the image matrix.
By keeping only the top k singular values and corresponding vectors, we can approximate the image.
This low-rank approximation preserves most visual details while using significantly fewer data values — achieving image compression.

Implementation Details
I have chosen the Full C Implementation (Option B) for this project.
All operations — image reading, SVD computation, reconstruction, and writing output — are performed in C without using Python.

Reading Image
The input image is read directly from a pgm file or jpg file.
I implemented a function to extract image dimensions and pixel values into an array.

Matrix Operations in C:
Wrote functions for basic matrix operations such as multiplication, transpose, normalization and for calculating forberian norm.

SVD Algorithm Implementation:
Implemented SVD using the Power Iteration Method to find eigenvalues and eigenvectors.
The algorithm extracts singular values by computing eigenvalues of ATA and then obtains the corresponding singular vectors.

Truncated SVD and Reconstruction:
Selected only the top k singular values.
Reconstructed the image matrix 

Writing Reconstructed Images:
The reconstructed matrix was converted back into pixel intensities (0–255) and saved as a new pgm file which then can be converted to jpg file.

Error Calculation:
Computed Frobenius error to measure how close Ak is to A
