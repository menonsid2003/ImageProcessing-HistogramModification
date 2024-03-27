#include <algorithm>
#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <string>
#include <sstream>
#include <cmath>
#include <cstdio>
#include <cassert>

using namespace std;

struct GeneralParameters
{
    string inputFile;
    string outputFile;
    // int numRois;
};

struct CommandParameters
{
    string command;
    vector<string> parameters;
};

struct PPM
{
    vector<vector<int>> Red;
    vector<vector<int>> Green;
    vector<vector<int>> Blue;
    int width = 0;
    int height = 0;
};

// Function to read parameters from a file
pair<GeneralParameters, vector<CommandParameters>> readParameters(const string &filename)
{
    ifstream file(filename);
    if (!file.is_open())
    {
        cerr << "Error opening file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    // Read the first line for general parameters
    GeneralParameters generalParams;
    file >> generalParams.inputFile >> generalParams.outputFile; // >> generalParams.numRois;

    // Read the rest of the lines for command parameters
    vector<CommandParameters> parameters;

    string line;
    while (getline(file, line))
    {
        if (line.empty())
            continue; // Ignore empty lines

        istringstream iss(line);
        CommandParameters commandParams;

        // Read the command
        iss >> commandParams.command;

        // Read parameters
        string param;
        while (iss >> param)
        {
            commandParams.parameters.push_back(param);
        }

        parameters.push_back(commandParams);
    }

    file.close();

    return {generalParams, parameters};
}

// Function to generate a histogram for a given ROI
vector<int> generateHistogram(const vector<vector<int>> &roi)
{
    vector<int> histogram(256, 0);

    for (const auto &row : roi)
    {
        for (int pixel : row)
        {
            histogram[pixel]++;
        }
    }

    return histogram;
}

pair<int, int> findMinMaxIndices(const vector<int> &array)
{
    int minIndex = -1; // Initialize to an invalid index
    int maxIndex = -1; // Initialize to an invalid index

    for (std::size_t i = 0; i < array.size(); ++i)
    {
        if (array[i] >= 1)
        {
            if (minIndex == -1)
            {
                // First occurrence of an element >= 1
                minIndex = i;
            }
            // Update maxIndex on each iteration to get the last occurrence
            maxIndex = i;
        }
    }

    return {minIndex, maxIndex};
}

// Function to print a histogram to the console
void printHistogram(const vector<int> &histogram)
{
    for (std::size_t i = 0; i < histogram.size(); ++i)
    {
        std::cout << i << ": " << histogram[i] << endl;
    }
}

// read from pgm (greyscale)
vector<vector<int>> readPGM(const string &filename)
{
    ifstream file(filename, ios::binary);

    if (!file.is_open())
    {
        cerr << "Error opening file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    string magicNumber; // magicnumber is the P2 or P5 at the top of a pgm file signifying whether its in ascii or binary.
    int width, height, maxPixelValue;

    // Read PGM header
    file >> magicNumber >> width >> height >> maxPixelValue;

    if (magicNumber != "P5")
    {
        cerr << "Invalid PGM file format: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    // Read pixel values into a 2D array
    vector<vector<int>> pixels(height, vector<int>(width));

    file.ignore(); // Ignore the newline character before pixel data

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            char pixelValue;
            file.read(&pixelValue, sizeof(char));

            // Convert binary character to unsigned int
            pixels[i][j] = bitset<8>(pixelValue).to_ulong();
        }
    }

    file.close();

    return pixels;
}

// write to pgm (greyscale)
void writePGM(const string &filename, const vector<vector<int>> &pixels)
{
    ofstream file(filename, ios::binary);

    if (!file.is_open())
    {
        cerr << "Error opening file for writing: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    const string magicNumber = "P5";
    const int width = pixels[0].size();
    const int height = pixels.size();
    const int maxPixelValue = 255; // limit pixel brightness to 255

    // Write PGM header
    file << magicNumber << "\n"
         << width << " " << height << "\n"
         << maxPixelValue << "\n";

    // Write pixel values
    for (const auto &row : pixels)
    {
        for (int value : row)
        {
            char pixelValue = static_cast<char>(value);
            file.write(&pixelValue, sizeof(char));
        }
    }

    file.close();
}

// flip image over the y axis
vector<vector<int>> flipImage(const vector<vector<int>> &inputImage)
{
    int height = inputImage.size();
    int width = inputImage[0].size();

    vector<vector<int>> flippedImage(height, vector<int>(width));

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            flippedImage[i][j] = inputImage[i][width - 1 - j];
        }
    }

    return flippedImage;
}

// Function to rotate an image either 90 or -90 degrees
vector<vector<int>> rotateImage(const vector<vector<int>> &inputImage, int degree)
{
    int height = inputImage.size();
    int width = inputImage[0].size();

    vector<vector<int>> rotatedImage;

    if (degree == 90 || degree == -90)
    {
        rotatedImage.resize(width, vector<int>(height));

        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                if (degree == 90)
                {
                    rotatedImage[j][height - 1 - i] = inputImage[i][j];
                }
                else // degree == -90
                {
                    rotatedImage[width - 1 - j][i] = inputImage[i][j];
                }
            }
        }
    }
    else if (degree == 180 || degree == -180)
    {
        rotatedImage.resize(height, vector<int>(width));

        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                rotatedImage[height - 1 - i][width - 1 - j] = inputImage[i][j];
            }
        }
    }
    else if (degree == 270 || degree == -270)
    {
        rotatedImage.resize(width, vector<int>(height));

        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                if (degree == 270)
                {
                    rotatedImage[width - 1 - j][i] = inputImage[i][j];
                }
                else // degree == -270
                {
                    rotatedImage[j][height - 1 - i] = inputImage[i][j];
                }
            }
        }
    }

    return rotatedImage;
}

PPM rotateImage(const PPM &inputImage, int degree)
{
    int height = inputImage.height;
    int width = inputImage.width;

    PPM rotatedImage;
    rotatedImage.width = (degree % 180 == 0) ? width : height;
    rotatedImage.height = (degree % 180 == 0) ? height : width;

    rotatedImage.Red.resize(rotatedImage.height, vector<int>(rotatedImage.width, 0));
    rotatedImage.Green.resize(rotatedImage.height, vector<int>(rotatedImage.width, 0));
    rotatedImage.Blue.resize(rotatedImage.height, vector<int>(rotatedImage.width, 0));

    if (degree == 90 || degree == -90)
    {
        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                if (degree == 90)
                {
                    rotatedImage.Red[j][height - 1 - i] = inputImage.Red[i][j];
                    rotatedImage.Green[j][height - 1 - i] = inputImage.Green[i][j];
                    rotatedImage.Blue[j][height - 1 - i] = inputImage.Blue[i][j];
                }
                else // degree == -90
                {
                    rotatedImage.Red[width - 1 - j][i] = inputImage.Red[i][j];
                    rotatedImage.Green[width - 1 - j][i] = inputImage.Green[i][j];
                    rotatedImage.Blue[width - 1 - j][i] = inputImage.Blue[i][j];
                }
            }
        }
    }
    else if (degree == 180 || degree == -180)
    {
        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                rotatedImage.Red[height - 1 - i][width - 1 - j] = inputImage.Red[i][j];
                rotatedImage.Green[height - 1 - i][width - 1 - j] = inputImage.Green[i][j];
                rotatedImage.Blue[height - 1 - i][width - 1 - j] = inputImage.Blue[i][j];
            }
        }
    }
    else if (degree == 270 || degree == -270)
    {
        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                if (degree == 270)
                {
                    rotatedImage.Red[width - 1 - j][i] = inputImage.Red[i][j];
                    rotatedImage.Green[width - 1 - j][i] = inputImage.Green[i][j];
                    rotatedImage.Blue[width - 1 - j][i] = inputImage.Blue[i][j];
                }
                else // degree == -270
                {
                    rotatedImage.Red[j][height - 1 - i] = inputImage.Red[i][j];
                    rotatedImage.Green[j][height - 1 - i] = inputImage.Green[i][j];
                    rotatedImage.Blue[j][height - 1 - i] = inputImage.Blue[i][j];
                }
            }
        }
    }

    return rotatedImage;
}

// Function to scale up an image without changing the image size, practically a zoom in
vector<vector<int>> scaleUpImage(const vector<vector<int>> &inputImage, double ratio)
{
    int size = inputImage.size();

    vector<vector<int>> scaledImage(size, vector<int>(size));

    for (int i = 0; i < size; ++i)
    {
        for (int j = 0; j < size; ++j)
        {
            int i2 = (int)floor((float)i / ratio);
            int j2 = (int)floor((float)j / ratio);
            scaledImage[i][j] = inputImage[i2][j2];
        }
    }

    return scaledImage;
}

// Function to subtract value from a pixel if it crosses the threshold
vector<vector<int>> subtractROI(const vector<vector<int>> &roi, int threshold, int value)
{
    vector<vector<int>> result = roi;

    for (auto &row : result)
    {
        for (int &pixel : row)
        {
            if (pixel < threshold)
            {
                pixel = max(0, pixel - value);
            }
        }
    }

    return result;
}

// Function to add a value to every pixel in an ROI
vector<vector<int>> addBrightnessToROI(const vector<vector<int>> &roi, int brightnessValue)
{
    vector<vector<int>> result = roi;

    for (auto &row : result)
    {
        for (int &pixel : row)
        {
            pixel = min(255, pixel + brightnessValue);
        }
    }

    return result;
}

// Function to take a given threshold value and set every pixel to black or white depending on if > or < than threshold
vector<vector<int>> binarizeROI(const vector<vector<int>> &roi, int threshold)
{
    vector<vector<int>> result = roi;

    for (auto &row : result)
    {
        for (int &pixel : row)
        {
            pixel = (pixel >= threshold) ? 255 : 0;
        }
    }

    return result;
}

// Function to select a square Region of Interest (ROI) from an image, clamped to image borders
vector<vector<int>> selectROI(const vector<vector<int>> &inputImage, int startX, int startY, int size)
{
    int imageHeight = inputImage.size();
    int imageWidth = inputImage[0].size();

    // Clamp startX and startY to be within valid range
    startX = max(0, min(startX, imageWidth - 1));
    startY = max(0, min(startY, imageHeight - 1));

    // Calculate the actual width and height of the ROI to ensure it stays within image boundaries
    int actualWidth = min(size, imageWidth - startX);
    int actualHeight = min(size, imageHeight - startY);

    vector<vector<int>> roiImage(actualHeight, vector<int>(actualWidth));

    for (int i = 0; i < actualHeight; ++i)
    {
        for (int j = 0; j < actualWidth; ++j)
        {
            roiImage[i][j] = inputImage[startY + i][startX + j];
        }
    }

    return roiImage;
}

// Specifically for localized histogram stretching to allow for width heigh input
vector<vector<int>> selectROI(const vector<vector<int>> &inputImage, int startX, int startY, int width, int height)
{
    int imageHeight = inputImage.size();
    int imageWidth = inputImage[0].size();

    // Clamp startX and startY to be within valid range
    startX = max(0, min(startX, imageWidth - 1));
    startY = max(0, min(startY, imageHeight - 1));

    // Calculate the actual width and height of the ROI to ensure it stays within image boundaries
    int actualWidth = min(width, imageWidth - startX);
    int actualHeight = min(height, imageHeight - startY);

    vector<vector<int>> roiImage(actualHeight, vector<int>(actualWidth));

    for (int i = 0; i < actualHeight; ++i)
    {
        for (int j = 0; j < actualWidth; ++j)
        {
            roiImage[i][j] = inputImage[startY + i][startX + j];
        }
    }

    return roiImage;
}

void replaceROI(vector<vector<int>> &originalImage, const vector<vector<int>> &roi, int startX, int startY)
{
    int roiHeight = roi.size();
    int roiWidth = roi[0].size();

    int imageHeight = originalImage.size();
    int imageWidth = originalImage[0].size();

    // Check if the given ROI coordinates are valid
    if (startX < 0 || startY < 0 || startX + roiWidth > imageWidth || startY + roiHeight > imageHeight)
    {
        cerr << "Invalid ROI coordinates" << endl;
        return;
    }

    // Replace the original elements with ROI elements
    for (int i = 0; i < roiHeight; ++i)
    {
        for (int j = 0; j < roiWidth; ++j)
        {
            originalImage[startY + i][startX + j] = roi[i][j];
        }
    }
}

/*                                                                      COLOR                                                                       */
// read from ppm (color)
PPM readPPM(const string &filename)
{
    std::fstream file(filename, fstream::in | fstream::binary);

    if (!file.is_open())
    {
        cerr << "Error opening file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    string magicNumber;
    int width, height, maxPixelValue;

    // Read PPM header
    file >> magicNumber >> width >> height >> maxPixelValue;

    if (magicNumber != "P6")
    {
        cerr << "Invalid PPM file format: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    // Ensure the maximum pixel value is within the expected range (0-255)
    if (maxPixelValue > 255)
    {
        cerr << "Unsupported max pixel value: " << maxPixelValue << endl;
        exit(EXIT_FAILURE);
    }

    PPM pixels;
    pixels.width = width;   // Assign width
    pixels.height = height; // Assign height

    // Initialize the Red, Green, and Blue vectors
    pixels.Red.resize(height, vector<int>(width, 0));
    pixels.Green.resize(height, vector<int>(width, 0));
    pixels.Blue.resize(height, vector<int>(width, 0));

    file.ignore(); // Ignore the newline character before pixel data

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                char pixelValue;
                file.read(&pixelValue, sizeof(char));

                // Convert binary character to unsigned int
                if (k == 0)
                {
                    pixels.Red[i][j] = bitset<8>(pixelValue).to_ulong();
                    // std::cout << bitset<8>(pixelValue).to_ulong() << " ";

                    // pixels.Red[i][j] = file.get();
                    //  cout << pixels.Red[i][j] << endl;
                }
                else if (k == 1)
                {
                    pixels.Green[i][j] = bitset<8>(pixelValue).to_ulong();
                    // pixels.Green[i][j] = file.get();
                    //  cout << pixels.Green[i][j] << endl;
                }
                else if (k == 2)
                {
                    pixels.Blue[i][j] = bitset<8>(pixelValue).to_ulong();
                    // pixels.Blue[i][j] = file.get();
                    //  cout << pixels.Blue[i][j] << endl;
                }
            }
        }
    }

    // vector<int> hist = generateHistogram(pixels.Red);
    // printHistogram(hist);
    // hist = generateHistogram(pixels.Green);
    // printHistogram(hist);
    // hist = generateHistogram(pixels.Blue);
    // printHistogram(hist);

    file.close();

    return pixels;
}

// write to ppm (color)
void writePPM(const string &filename, const PPM &pixels)
{
    ofstream file(filename, ios::binary);

    if (!file.is_open())
    {
        cerr << "Error opening file for writing: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    const string magicNumber = "P6";
    const int width = pixels.width;
    const int height = pixels.height;
    const int maxPixelValue = 255; // limit brightness to 255

    // Write PPM header
    file << magicNumber << "\n"
         << width << " " << height << "\n"
         << maxPixelValue << "\n";

    // Write pixel values (3 for loops to account for color values)
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            // cout << pixels.Red[i][j] << endl;
            //  Save each channel into char and write into image
            char pixelValue = static_cast<char>(pixels.Red[i][j]);
            file.write(&pixelValue, sizeof(char));
            pixelValue = static_cast<char>(pixels.Green[i][j]);
            file.write(&pixelValue, sizeof(char));
            pixelValue = static_cast<char>(pixels.Blue[i][j]);
            file.write(&pixelValue, sizeof(char));
        }
    }

    file.close();
}

// select roi from ppm
PPM selectROI(const PPM &inputImage, int startX, int startY, int size)
{
    int imageHeight = inputImage.height;
    int imageWidth = inputImage.width;

    // Clamp startX and startY to be within valid range
    startX = max(0, min(startX, imageWidth - 1));
    startY = max(0, min(startY, imageHeight - 1));

    // Calculate the actual width and height of the ROI to ensure it stays within image boundaries
    int actualWidth = min(size, imageWidth - startX);
    int actualHeight = min(size, imageHeight - startY);

    PPM roiImage;
    roiImage.width = actualWidth;
    roiImage.height = actualHeight;
    roiImage.Red.resize(actualHeight, vector<int>(actualWidth));
    roiImage.Green.resize(actualHeight, vector<int>(actualWidth));
    roiImage.Blue.resize(actualHeight, vector<int>(actualWidth));

    for (int i = 0; i < actualHeight; ++i)
    {
        for (int j = 0; j < actualWidth; ++j)
        {
            roiImage.Red[i][j] = inputImage.Red[startY + i][startX + j];
            roiImage.Green[i][j] = inputImage.Green[startY + i][startX + j];
            roiImage.Blue[i][j] = inputImage.Blue[startY + i][startX + j];
        }
    }

    return roiImage;
}

// flip ppm
vector<vector<vector<int>>> flipImage(const vector<vector<vector<int>>> &inputImage)
{
    int height = inputImage.size();
    int width = inputImage[0].size();

    vector<vector<vector<int>>> flippedImage(height, vector<vector<int>>(width, vector<int>(3)));

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            flippedImage[i][j] = inputImage[i][width - 1 - j];
        }
    }

    return flippedImage;
}

// rotate ppm
vector<vector<vector<int>>> rotateImage(const vector<vector<vector<int>>> &inputImage, int degree)
{
    int height = inputImage.size();
    int width = inputImage[0].size();

    vector<vector<vector<int>>> rotatedImage;

    if (degree == 90 || degree == -90)
    {
        rotatedImage.resize(width, vector<vector<int>>(height, vector<int>(3)));

        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                if (degree == 90)
                {
                    rotatedImage[j][height - 1 - i] = inputImage[i][j];
                }
                else // degree == -90
                {
                    rotatedImage[width - 1 - j][i] = inputImage[i][j];
                }
            }
        }
    }
    else if (degree == 180 || degree == -180)
    {
        rotatedImage.resize(height, vector<vector<int>>(width, vector<int>(3)));

        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                rotatedImage[height - 1 - i][width - 1 - j] = inputImage[i][j];
            }
        }
    }
    else if (degree == 270 || degree == -270)
    {
        rotatedImage.resize(width, vector<vector<int>>(height, vector<int>(3)));

        for (int i = 0; i < height; ++i)
        {
            for (int j = 0; j < width; ++j)
            {
                if (degree == 270)
                {
                    rotatedImage[width - 1 - j][i] = inputImage[i][j];
                }
                else // degree == -270
                {
                    rotatedImage[j][height - 1 - i] = inputImage[i][j];
                }
            }
        }
    }

    return rotatedImage;
}

// multiplicative color mod
vector<vector<vector<int>>> multiplicativeColor(const vector<vector<vector<int>>> &roi, double moreC)
{
    int size = roi.size();

    vector<vector<vector<int>>> multRoi(size, vector<vector<int>>(size, vector<int>(3)));

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            multRoi[i][j][0] = min(255, max(0, roi[i][j][0] * static_cast<int>(moreC)));
            multRoi[i][j][1] = min(255, max(0, roi[i][j][1] * static_cast<int>(moreC)));
            multRoi[i][j][2] = min(255, max(0, roi[i][j][2] * static_cast<int>(moreC)));
        }
    }

    return multRoi;
}

// additive color mod with clamping
vector<vector<vector<int>>> additiveColor(const vector<vector<vector<int>>> &roi, double AC)
{
    int size = roi.size();

    vector<vector<vector<int>>> result(size, vector<vector<int>>(size, vector<int>(3)));

    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            // Add the brightness with clamping
            result[i][j][0] = min(255, max(0, roi[i][j][0] + static_cast<int>(AC)));
            result[i][j][1] = min(255, max(0, roi[i][j][1] + static_cast<int>(AC)));
            result[i][j][2] = min(255, max(0, roi[i][j][2] + static_cast<int>(AC)));
        }
    }

    return result;
}

void writeRedPPM(const string &filename, int width, int height)
{
    ofstream file(filename, ios::binary);

    if (!file.is_open())
    {
        cerr << "Error opening file for writing: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    const string magicNumber = "P6";
    const int maxPixelValue = 255; // limit brightness to 255

    // Write PPM header
    file << magicNumber << "\n"
         << width << " " << height << "\n"
         << maxPixelValue << "\n";

    // Write pure red pixel values
    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            char redPixel = static_cast<char>(100); // Red channel at full brightness
            char greenPixel = static_cast<char>(0); // Green channel at zero brightness
            char bluePixel = static_cast<char>(0);  // Blue channel at zero brightness

            // Write each color channel separately
            file.write(&redPixel, sizeof(char));
            file.write(&greenPixel, sizeof(char));
            file.write(&bluePixel, sizeof(char));
        }
    }

    file.close();
}

// Function to perform histogram stretching on a grey-level image
vector<vector<int>> stretchHistogram(const vector<vector<int>> &inputImage, int A, int B)
{
    int height = inputImage.size();
    int width = (height > 0) ? inputImage[0].size() : 0;

    if (width == 0)
    {
        // Handle empty inputImage
        return {};
    }

    vector<int> hist = generateHistogram(inputImage);
    // printHistogram(hist);
    pair<int, int> minMax = findMinMaxIndices(hist);
    int I_min = minMax.first;
    int I_max = minMax.second;
    double c = 1.05 * I_min;
    double d = 0.95 * I_max;

    if (c == d)
    {
        // std::cerr << "C = D error" << endl;
        //  std::terminate();
    }

    vector<vector<int>> stretchedImage(height, vector<int>(width, 0));

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            int stretchedValue;

            if (inputImage[i][j] < c)
            {
                stretchedValue = A;
            }
            else if (inputImage[i][j] > d)
            {
                stretchedValue = B;
            }
            else
            {
                stretchedValue = static_cast<int>((inputImage[i][j] - c) * ((B - A) / (d - c)) + A);
            }

            stretchedImage[i][j] = min(255, max(0, stretchedValue));
        }
    }

    return stretchedImage;
}

// Function to apply histogram stretching to a local window around a pixel
vector<vector<int>> slidingWindowStretch(const vector<vector<int>> &inputImage, int A, int B, int S)
{
    int height = inputImage.size();
    int width = inputImage[0].size();

    vector<vector<int>> stretchedImage(height, vector<int>(width));

    for (int i = 0; i < height; ++i)
    {
        for (int j = 0; j < width; ++j)
        {
            // Define the boundaries of the local window
            int startY = max(0, i - S / 2);
            int endY = min(height - 1, i + S / 2);
            int startX = max(0, j - S / 2);
            int endX = min(width - 1, j + S / 2);

            // Extract the local window
            vector<vector<int>> localWindow;
            for (int y = startY; y <= endY; ++y)
            {
                vector<int> row;
                for (int x = startX; x <= endX; ++x)
                {
                    row.push_back(inputImage[y][x]);
                }
                localWindow.push_back(row);
            }

            // Apply histogram stretching to the local window
            vector<vector<int>> stretchedWindow = stretchHistogram(localWindow, A, B);

            // Update the pixel value in the result image
            stretchedImage[i][j] = stretchedWindow[S / 2][S / 2]; // Use the center pixel of the local window
        }
    }

    return stretchedImage;
}

vector<vector<int>> localizedQuadrantStretch(const vector<vector<int>> &inputImage, int A, int B)
{
    int height = inputImage.size();
    int width = inputImage[0].size();

    // Split the input image into four quadrants
    int midRow = height / 2;
    int midCol = width / 2;

    vector<vector<int>> upperLeftQuadrant = selectROI(inputImage, 0, 0, midCol, midRow);
    vector<vector<int>> upperRightQuadrant = selectROI(inputImage, midCol, 0, width - midCol, midRow);
    vector<vector<int>> lowerLeftQuadrant = selectROI(inputImage, 0, midRow, midCol, height - midRow);
    vector<vector<int>> lowerRightQuadrant = selectROI(inputImage, midCol, midRow, width - midCol, height - midRow);

    // Apply histogram stretching to each quadrant
    vector<vector<int>> stretchedUpperLeft = stretchHistogram(upperLeftQuadrant, A, B);
    vector<vector<int>> stretchedUpperRight = stretchHistogram(upperRightQuadrant, A, B);
    vector<vector<int>> stretchedLowerLeft = stretchHistogram(lowerLeftQuadrant, A, B);
    vector<vector<int>> stretchedLowerRight = stretchHistogram(lowerRightQuadrant, A, B);

    // Combine the results into the final stretched image
    vector<vector<int>> stretchedImage(height, vector<int>(width));

    for (int i = 0; i < midRow; ++i)
    {
        for (int j = 0; j < midCol; ++j)
        {
            stretchedImage[i][j] = stretchedUpperLeft[i][j];
            stretchedImage[i][j + midCol] = stretchedUpperRight[i][j];
            stretchedImage[i + midRow][j] = stretchedLowerLeft[i][j];
            stretchedImage[i + midRow][j + midCol] = stretchedLowerRight[i][j];
        }
    }

    return stretchedImage;
}

PPM colorHistogramStretch(const PPM &inputImage, char channel, int A, int B)
{
    PPM stretchedImage;
    stretchedImage.height = inputImage.height;
    stretchedImage.width = inputImage.width;
    stretchedImage.Red.resize(inputImage.width, vector<int>(inputImage.height));
    stretchedImage.Red = inputImage.Red;
    stretchedImage.Green.resize(inputImage.width, vector<int>(inputImage.height));
    stretchedImage.Green = inputImage.Green;
    stretchedImage.Blue.resize(inputImage.width, vector<int>(inputImage.height));
    stretchedImage.Blue = inputImage.Blue;

    // If channel is 'A', apply stretching to all channels
    if (channel == 'A' || channel == 'a')
    {
        stretchedImage.Red = localizedQuadrantStretch(inputImage.Red, A, B);
        stretchedImage.Green = localizedQuadrantStretch(inputImage.Green, A, B);
        stretchedImage.Blue = localizedQuadrantStretch(inputImage.Blue, A, B);
    }
    // If channel is 'R', apply stretching to the Red channel only
    else if (channel == 'R' || channel == 'r')
    {
        stretchedImage.Red = localizedQuadrantStretch(inputImage.Red, A, B);
    }
    // If channel is 'G', apply stretching to the Green channel only
    else if (channel == 'G' || channel == 'g')
    {
        stretchedImage.Green = localizedQuadrantStretch(inputImage.Green, A, B);
    }
    // If channel is 'B', apply stretching to the Blue channel only
    else if (channel == 'B' || channel == 'b')
    {
        stretchedImage.Blue = localizedQuadrantStretch(inputImage.Blue, A, B);
    }
    // If an invalid channel is provided, return the original image
    else
    {
        cerr << "Invalid channel. No stretching applied." << endl;
        return inputImage;
    }

    return stretchedImage;
}

void performColorStretch(const std::string &colorFile, char channel, int startX, int startY, int size, int A, int B)
{
    static int roiColorNum = 0; // Static variable to keep track of roiColorNum

    // Read PPM image
    PPM ppmImage = readPPM(colorFile);

    // Select ROI from the PPM image
    PPM roiColor = selectROI(ppmImage, startX, startY, size);

    // Save the original ROI color
    std::string originalPath = "histogram/color/original/";
    writePPM(originalPath + "Roi" + to_string(roiColorNum) + "_original.ppm", roiColor);

    // Set folder path based on the channel
    std::string folderpath;
    switch (channel)
    {
    case 'R':
        folderpath = "histogram/color/red/";
        break;
    case 'B':
        folderpath = "histogram/color/blue/";
        break;
    case 'G':
        folderpath = "histogram/color/green/";
        break;
    case 'A':
        folderpath = "histogram/color/combined/";
        break;
    default:
        std::cerr << "Invalid channel specified." << std::endl;
        return;
    }

    // Print start message
    std::cout << channel << " Start" << std::endl;

    // Perform color histogram stretching
    PPM roiStretch = colorHistogramStretch(roiColor, channel, A, B);

    // Print end message
    std::cout << channel << " End" << std::endl;

    // Write the stretched image to a file
    writePPM(folderpath + "Roi" + to_string(roiColorNum) + "_stretched.ppm", roiStretch);

    // Rotate the images
    PPM modded;

    modded = rotateImage(roiColor, 90);
    writePPM(originalPath + "Roi" + to_string(roiColorNum) + "_original_rotate90.ppm", modded);

    modded = rotateImage(roiColor, 180);
    writePPM(originalPath + "Roi" + to_string(roiColorNum) + "_original_rotate180.ppm", modded);

    modded = rotateImage(roiColor, 270);
    writePPM(originalPath + "Roi" + to_string(roiColorNum) + "_original_rotate270.ppm", modded);

    modded = rotateImage(roiStretch, 90);
    writePPM(folderpath + "Roi" + to_string(roiColorNum) + "_stretched_rotate90.ppm", modded);

    modded = rotateImage(roiStretch, 180);
    writePPM(folderpath + "Roi" + to_string(roiColorNum) + "_stretched_rotate180.ppm", modded);

    modded = rotateImage(roiStretch, 270);
    writePPM(folderpath + "Roi" + to_string(roiColorNum) + "_stretched_rotate270.ppm", modded);

    // Print completion message
    std::cout << "Stretched image saved as Roi" << to_string(roiColorNum) << "stretched.ppm in " << folderpath << std::endl;
    // Increment roiColorNum for the next iteration
    roiColorNum++;
}

int main()
{
    string parametersFilename = "parameters.txt";
    pair<GeneralParameters, vector<CommandParameters>> paramsPair = readParameters(parametersFilename);

    // Process general parameters
    GeneralParameters generalParams = paramsPair.first;
    std::cout << "Input File: " << generalParams.inputFile << endl;
    std::cout << "Output File: " << generalParams.outputFile << endl;
    // std::cout << "Number of ROIs: " << generalParams.numRois << endl;

    vector<vector<int>> image = readPGM(generalParams.inputFile);

    vector<vector<vector<int>>> allROIs;
    int roinum = 0;
    int roiColorNum = 0;

    //   iterate through each line and its parameters
    vector<CommandParameters> commandParams = paramsPair.second;
    int x = 0;

    for (const auto &cmdParams : commandParams)
    {
        std::cout << "Command21: " << cmdParams.command << endl;

        if (cmdParams.command == "//")
        {
            cout << "Line Commented out, skipping.." << endl;
            continue;
        }
        else if (cmdParams.command == "sub")
        {
            if (cmdParams.parameters.size() == 5)
            {
                int startX = stoi(cmdParams.parameters[0]);
                int startY = stoi(cmdParams.parameters[1]);
                int size = stoi(cmdParams.parameters[2]);
                int threshold = stoi(cmdParams.parameters[3]);
                int value = stoi(cmdParams.parameters[4]);

                vector<vector<int>> Roi1 = selectROI(image, startX, startY, size);
                allROIs.push_back(Roi1);

                vector<vector<int>> deepCopyRoi1 = Roi1;
                deepCopyRoi1 = subtractROI(deepCopyRoi1, threshold, value);
                replaceROI(image, deepCopyRoi1, startX, startY);
                writePGM(generalParams.outputFile, image);
            }
            else
            {
                cerr << "Invalid parameters for command: " << cmdParams.command << endl;
            }
        }
        else if (cmdParams.command == "add")
        {
            if (cmdParams.parameters.size() == 4)
            {
                int startX = stoi(cmdParams.parameters[0]);
                int startY = stoi(cmdParams.parameters[1]);
                int size = stoi(cmdParams.parameters[2]);
                int value = stoi(cmdParams.parameters[3]);

                vector<vector<int>> Roi2 = selectROI(image, startX, startY, size);
                allROIs.push_back(Roi2);

                vector<vector<int>> deepCopyRoi2 = Roi2;
                deepCopyRoi2 = addBrightnessToROI(deepCopyRoi2, value);
                replaceROI(image, deepCopyRoi2, startX, startY);
                writePGM(generalParams.outputFile, image);
            }
            else
            {
                cerr << "Invalid parameters for command: " << cmdParams.command << endl;
            }
        }
        else if (cmdParams.command == "binarize")
        {
            if (cmdParams.parameters.size() == 4)
            {
                int startX = stoi(cmdParams.parameters[0]);
                int startY = stoi(cmdParams.parameters[1]);
                int size = stoi(cmdParams.parameters[2]);
                int threshold = stoi(cmdParams.parameters[3]);

                vector<vector<int>> Roi3 = selectROI(image, startX, startY, size);
                allROIs.push_back(Roi3);

                vector<vector<int>> deepCopyRoi3 = Roi3;
                deepCopyRoi3 = binarizeROI(deepCopyRoi3, threshold);
                replaceROI(image, deepCopyRoi3, startX, startY);
                writePGM(generalParams.outputFile, image);
            }
            else
            {
                cerr << "Invalid parameters for command: " << cmdParams.command << endl;
            }
        }
        else if (cmdParams.command == "*") //* indicates start of task 2 and following parameters such as scale ratio and brightness value
        {

            int i = 1;
            for (const auto &currRoi : allROIs)
            {
                // b1
                vector<vector<int>> modded = flipImage(currRoi);
                string folderpath = "Roi" + to_string(i) + "/" + "flipped/";
                writePGM(folderpath + "Roi" + to_string(i) + "flipped" + ".pgm", modded);
                modded = rotateImage(currRoi, 90);
                writePGM(folderpath + "Roi" + to_string(i) + "rotate90" + ".pgm", modded);
                modded = rotateImage(currRoi, -90);
                writePGM(folderpath + "Roi" + to_string(i) + "rotate-90" + ".pgm", modded);

                // b2
                folderpath = "Roi" + to_string(i) + "/" + "scaled/";
                vector<vector<int>> moddedScaled = currRoi;
                moddedScaled = scaleUpImage(currRoi, stoi(cmdParams.parameters[0]));
                writePGM(folderpath + "Roi" + to_string(i) + "Scaled" + ".pgm", moddedScaled);
                modded = flipImage(moddedScaled);
                writePGM(folderpath + "Roi" + to_string(i) + "Scaled" + "flipped" + ".pgm", modded);
                modded = rotateImage(moddedScaled, 90);
                writePGM(folderpath + "Roi" + to_string(i) + "Scaled" + "rotate90" + ".pgm", modded);
                modded = rotateImage(moddedScaled, -90);
                writePGM(folderpath + "Roi" + to_string(i) + "Scaled" + "rotate-90" + ".pgm", modded);

                // b3
                folderpath = "Roi" + to_string(i) + "/" + "brightness/";
                vector<vector<int>> moddedBright = currRoi;
                moddedBright = addBrightnessToROI(currRoi, stoi(cmdParams.parameters[1]));
                writePGM(folderpath + "Roi" + to_string(i) + "Brightened" + ".pgm", moddedBright);
                modded = flipImage(moddedBright);
                writePGM(folderpath + "Roi" + to_string(i) + "Brightened" + "flipped" + ".pgm", modded);
                modded = rotateImage(moddedBright, 90);
                writePGM(folderpath + "Roi" + to_string(i) + "Brightened" + "rotate90" + ".pgm", modded);
                modded = rotateImage(moddedBright, -90);
                writePGM(folderpath + "Roi" + to_string(i) + "Brightened" + "rotate-90" + ".pgm", modded);

                i++;
            }
        }
        else if (cmdParams.command.size() >= 4 && cmdParams.command.compare(cmdParams.command.size() - 4, 4, ".ppm") == 0)
        // a command ending in .ppm indicated a .ppm filename to be read in
        {

            if (cmdParams.parameters.size() == 5)
            {
                /*
                string ppmFilename = cmdParams.command;
                int startX = stoi(cmdParams.parameters[0]);
                int startY = stoi(cmdParams.parameters[1]);
                int size = stoi(cmdParams.parameters[2]);
                double multC = stod(cmdParams.parameters[3]);
                int addC = stoi(cmdParams.parameters[4]);
                // std::cout << multC;
                vector<vector<vector<int>>> ppmImage = readPPM(ppmFilename);

                vector<vector<vector<int>>> roi = selectROI(ppmImage, startX, startY, size);

                vector<vector<vector<int>>> multRoi = multiplicativeColor(roi, multC);
                vector<vector<vector<int>>> addRoi = additiveColor(roi, addC);

                vector<vector<vector<int>>> modded = multRoi;
                string folderpath = "color/multiplicative/";
                writePPM(folderpath + "multC.ppm", modded);
                modded = flipImage(multRoi);
                writePPM(folderpath + "multC" + "flipped" + ".ppm", modded);
                modded = rotateImage(multRoi, 90);
                writePPM(folderpath + "multC" + "rotate90" + ".ppm", modded);
                modded = rotateImage(multRoi, -90);
                writePPM(folderpath + "multC" + "rotate-90" + ".ppm", modded);

                folderpath = "color/additive/";
                writePPM(folderpath + "addC.ppm", addRoi);
                modded = flipImage(addRoi);
                writePPM(folderpath + "addC" + "flipped" + ".ppm", modded);
                modded = rotateImage(addRoi, 90);
                writePPM(folderpath + "addC" + "rotate90" + ".ppm", modded);
                modded = rotateImage(addRoi, -90);
                writePPM(folderpath + "addC" + "rotate-90" + ".ppm", modded);
                writePPM("addC.ppm", addRoi);
                */
            }
            else
            {
                cerr << "Invalid parameters for command: " << cmdParams.command << endl;
            }
        }
        else if (cmdParams.command == "allstretch")
        {
            if (cmdParams.parameters.size() == 5)
            {
                int startX = stoi(cmdParams.parameters[0]);
                int startY = stoi(cmdParams.parameters[1]);
                int size = stoi(cmdParams.parameters[2]);
                int A = static_cast<int>(stoi(cmdParams.parameters[3]));
                int B = static_cast<int>(stoi(cmdParams.parameters[4]));

                vector<vector<int>> roiStretch = selectROI(image, startX, startY, size);

                // Ca
                string folderpath = "histogram/original/";
                writePGM(folderpath + "Roi" + to_string(roinum) + "original" + ".pgm", roiStretch);

                vector<vector<int>> modded = rotateImage(roiStretch, 90);
                writePGM(folderpath + "Roi" + to_string(roinum) + "original" + "rotate90" + ".pgm", modded);

                modded = rotateImage(roiStretch, 180);
                writePGM(folderpath + "Roi" + to_string(roinum) + "original" + "rotate180" + ".pgm", modded);

                modded = rotateImage(roiStretch, 270);
                writePGM(folderpath + "Roi" + to_string(roinum) + "original" + "rotate270" + ".pgm", modded);

                // Cb
                folderpath = "histogram/stretch/";
                vector<vector<int>> stretchedRoi = stretchHistogram(roiStretch, A, B);
                writePGM(folderpath + "Roi" + to_string(roinum) + "Stretched" + ".pgm", stretchedRoi);

                modded = rotateImage(stretchedRoi, 90);
                writePGM(folderpath + "Roi" + to_string(roinum) + "Stretched" + "rotate90" + ".pgm", modded);

                modded = rotateImage(stretchedRoi, 180);
                writePGM(folderpath + "Roi" + to_string(roinum) + "Stretched" + "rotate180" + ".pgm", modded);

                modded = rotateImage(stretchedRoi, 270);
                writePGM(folderpath + "Roi" + to_string(roinum) + "Stretched" + "rotate270" + ".pgm", modded);

                // Cc
                folderpath = "histogram/localstretch/";
                vector<vector<int>> localStretch = roiStretch;

                localStretch = localizedQuadrantStretch(roiStretch, A, B);
                writePGM(folderpath + "Roi" + to_string(roinum) + "localStretch" + ".pgm", localStretch);

                modded = rotateImage(localStretch, 90);
                writePGM(folderpath + "Roi" + to_string(roinum) + "localStretch" + "rotate90" + ".pgm", modded);

                modded = rotateImage(localStretch, 180);
                writePGM(folderpath + "Roi" + to_string(roinum) + "localStretch" + "rotate180" + ".pgm", modded);

                modded = rotateImage(localStretch, 270);
                writePGM(folderpath + "Roi" + to_string(roinum) + "localStretch" + "rotate270" + ".pgm", modded);
                roinum++;
            }
        }
        else if (cmdParams.command == "slidingwindow")
        {
            if (cmdParams.parameters.size() == 6)
            {
                int startX = stoi(cmdParams.parameters[0]);
                int startY = stoi(cmdParams.parameters[1]);
                int size = stoi(cmdParams.parameters[2]);
                int A = static_cast<int>(stoi(cmdParams.parameters[3]));
                int B = static_cast<int>(stoi(cmdParams.parameters[4]));
                int S = static_cast<int>(stoi(cmdParams.parameters[5]));

                vector<vector<int>> roiStretch = selectROI(image, startX, startY, size);
                roiStretch = slidingWindowStretch(roiStretch, A, B, S);
                // Ca
                string folderpath = "histogram/";
                writePGM(folderpath + "Roi" + "SlidingWindow" + ".pgm", roiStretch);
            }
        }
        else if (cmdParams.command == "colorstretch")
        {
            if (cmdParams.parameters.size() == 7)
            {
                string colorFile = cmdParams.parameters[0];
                char channel = cmdParams.parameters[1][0];
                int startX = stoi(cmdParams.parameters[2]);
                int startY = stoi(cmdParams.parameters[3]);
                int size = stoi(cmdParams.parameters[4]);
                int A = static_cast<int>(stoi(cmdParams.parameters[5]));
                int B = static_cast<int>(stoi(cmdParams.parameters[6]));

                performColorStretch(colorFile, channel, startX, startY, size, A, B);
            }
            roiColorNum++;
        }
        else
        {
            cerr << "Unknown command: " << cmdParams.command << endl;
        }
        x++;
    }
    std::cout << "Exiting.." << endl;
    return 0;
}