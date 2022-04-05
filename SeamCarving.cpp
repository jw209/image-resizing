#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

// using declarations
using std::cout;
using std::endl;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::stringstream;
using std::string;
using std::size_t;
using std::stoi;
using std::atoi;
using std::abs;
using std::min;

// functions declarations
vector<vector<int>> transposeMatrix(vector<vector<int>> matrix);
vector<vector<int>> createEnergyMatrix(vector<vector<int>> mat);
int getEnergyOfIndex(int left, int right, int top, int bottom, int numAtIndex, int selection);
vector<vector<int>> createCumulativeEnergyMatrix(vector<vector<int>> energyMap);
vector<int> findSeam(vector<vector<int>> cumulativeEnergyMatrix);
void removeSeams(vector<vector<int>> &matrix, vector<int> seam);
vector<vector<int>> removeSeamsDriver(vector<vector<int>> mat, int numV, int numH);

int main (int argc, char* argv[]) {

    // initialize file input variables
    if (argv[1] && argv[2] && argv[3]) {

        /*
        // This section of code gets processes the .pgm file and outputs its information
        */
        int numV = atoi(argv[2]);
        int numH = atoi(argv[3]);

        ifstream infile(argv[1]);
        stringstream ss;
        string input = "";

        // get version of pgm file
        getline(infile, input);
        if(input.compare("P2") != 0) {
            cout << "Incorrect .pgm version" << endl;
        } else {
            cout << "Version: " << input << endl;
        }

        // get comment
        getline(infile, input);
        cout << "Comment: " << input << endl;

        // process rows and columns
        getline(infile, input);
        size_t spacePos = input.find(" ");
        size_t endPos = input.find("\n");
        string colCount = input.substr(0, spacePos);
        string rowCount = input.substr(spacePos+1, endPos);

        // row and column count to integers
        int rows = stoi(rowCount);
        int columns = stoi(colCount);
        cout << "Columns and rows (respectively): " << columns << " " << rows << endl;

        // get the scale of the greyscale
        getline(infile, input);
        int scale = stoi(input);
        cout << "Scale: " << scale << endl;

        /*
        // Read in greyscale map to vector
        */
        ss << infile.rdbuf();

        vector<vector<int>> mat(rows, vector<int>(columns));

        for (int row = 0; row < rows; ++row)
            for (int col = 0; col < columns; ++col)
                ss >> mat[row][col];

        infile.close();

        vector<vector<int>> outputMat(rows-numH, vector<int>(columns-numV));
        outputMat = removeSeamsDriver(mat, numV, numH);

        //UNCOMMENT TO OUTPUT FILE
        string fname = argv[1];
        size_t periodIndex = fname.find(".");
        string filename = fname.substr(0, periodIndex);

        string outfilename = filename+"_processedX_"+argv[2]+"_"+argv[3]+".pgm";
        ofstream outfile(outfilename);

        outfile << "P2\n";
        outfile << "# Processed by Jesse Wood's seam carving algorithm\n";
        outfile << outputMat[0].size() << " " << outputMat.size() << "\n";
        outfile << "255\n";

        for (int row = 0; row < outputMat.size(); ++row)
            for (int col = 0; col < outputMat[row].size(); ++col)
                outfile << outputMat[row][col] << " ";

        outfile.close();

    } else {

        // output if incorrect command is given
        cout << "Incorrect command. To run this program use \"./main 'filepath' 'vertical seams' 'horizontal seams'\"" << endl;
        return 1;
    }  

    return 0;
}

// TRANSPOSE MATRIX FUNCTION
// input: matrix to transpose
// output: transposed matrix
vector<vector<int>> transposeMatrix(vector<vector<int>> matrix) {

    std::vector<std::vector<int>> transposedMat(matrix[0].size(), std::vector<int>(matrix.size()));

    for (int row = 0; row < matrix.size(); ++row)
        for (int col = 0; col < matrix[row].size(); ++col)
            transposedMat[col][row] = matrix[row][col];
    
    return transposedMat;
}

// CREATE ENERGY MATRIX FUNCTION
// input: matrix from which to create an energy matrix
// output: energy matrix
vector<vector<int>> createEnergyMatrix(vector<vector<int>> mat) {

    vector<vector<int>> energyMap(mat.size(), vector<int>(mat[0].size()));

    int numAtIndex = 0, bottom = 0, right = 0, left = 0, top = 0;
    
    for (int i = 0; i < mat.size(); i++) {
        for (int j = 0; j < mat[i].size(); j++) {
            numAtIndex = mat[i][j];
            if (i == 0) {
                if (j == 0 || j == mat[i].size()-1) {
                    // process the middle of first row
                    //
                    // get values at relevant indices
                    if (j == 0) {
                        // process top left corner
                        bottom = mat[i+1][j];
                        right = mat[i][j+1];
                        energyMap[i][j] = getEnergyOfIndex(0, right, 0, bottom, numAtIndex, 1);
                    } else if (j == mat[i].size()-1) {
                        // process top right corner
                        left = mat[i][j-1];
                        bottom = mat[i+1][j];
                        energyMap[i][j] = getEnergyOfIndex(left, 0, 0, bottom, numAtIndex, 3);
                    } 
                } else {
                        left = mat[i][j-1];
                        right = mat[i][j+1];
                        bottom = mat[i+1][j];
                        energyMap[i][j] = getEnergyOfIndex(left, right, 0, bottom, numAtIndex, 2);
                }
            } else if (i == mat.size()-1) {
                if (j == 0 || j == mat[i].size()-1) {
                    if (j == 0) {
                        // process bottom left corner
                        right = mat[i][j+1];
                        top = mat[i-1][j];
                        energyMap[i][j] = getEnergyOfIndex(0, right, top, 0, numAtIndex, 7);
                    } else if (j == mat[i].size()-1) {
                        // process bottom right corner
                        top = mat[i-1][j];
                        left = mat[i][j-1];
                        energyMap[i][j] = getEnergyOfIndex(left, 0, top, 0, numAtIndex, 9);
                    } 
                } else {
                        // process the middle of last row
                        left = mat[i][j-1];
                        right = mat[i][j+1];
                        top = mat[i-1][j];
                        energyMap[i][j] = getEnergyOfIndex(left, right, top, 0, numAtIndex, 8);
                }
            } else {
                if (j == 0 || j == mat[i].size()-1) {
                    if (j == 0) {
                        // process farthest left part of matrix
                        right = mat[i][j+1];
                        top = mat[i-1][j];
                        bottom = mat[i+1][j];
                        energyMap[i][j] = getEnergyOfIndex(0, right, top, bottom, numAtIndex, 4);
                    } else if (j == mat[i].size()-1) {
                        // process farthest right part of matrix
                        left = mat[i][j-1];
                        top = mat[i-1][j];
                        bottom = mat[i+1][j];
                        energyMap[i][j] = getEnergyOfIndex(left, 0, top, bottom, numAtIndex, 6);
                    }
                } else {
                        // process inner parts of the matrix
                        left = mat[i][j-1];
                        right = mat[i][j+1];
                        top = mat[i-1][j];
                        bottom = mat[i+1][j];
                        energyMap[i][j] = getEnergyOfIndex(left, right, top, bottom, numAtIndex, 5);
                }
            }
        }
    }
    return energyMap;
}

// GET ENERGY OF INDEX FUNCTION
// input: value left of index, value right of index, value above index, value below index, value of index, selector
// output: energy at index
int getEnergyOfIndex(int left, int right, int top, int bottom, int numAtIndex, int selection) {

    int leftDiff = abs(numAtIndex - left);
    int rightDiff = abs(numAtIndex - right);
    int topDiff = abs(numAtIndex - top);
    int bottomDiff = abs(numAtIndex - bottom);

    switch (selection) {
        case 1: return (rightDiff + bottomDiff);
                break;
        case 2: return (leftDiff + rightDiff + bottomDiff);
                break;
        case 3: return (leftDiff + bottomDiff);
                break;
        case 4: return (topDiff + rightDiff + bottomDiff);
                break;
        case 5: return (topDiff + leftDiff + rightDiff + bottomDiff);
                break;
        case 6: return (topDiff + leftDiff + bottomDiff);
                break;
        case 7: return (topDiff + rightDiff);
                break;
        case 8: return (leftDiff + topDiff + rightDiff);
                break;
        case 9: return (leftDiff + topDiff);
                break;
        default: return 0;
    }
}

// CREATE CUMULATIVE ENERGY MATRIX FUNCTION
// input: energy map
// output: cumulative energy map
vector<vector<int>> createCumulativeEnergyMatrix(vector<vector<int>> energyMap) {

    vector<vector<int>> cumulativeEnergyMatrix(energyMap.size(), vector<int>(energyMap[0].size()));

    // load first last row with values from the energy map
    for (int j = 0; j < energyMap[0].size(); j++) {
        cumulativeEnergyMatrix[energyMap.size()-1][j] = energyMap[energyMap.size()-1][j];
    }

    // populate the rest of the matrix
    for (int i = energyMap.size()-2; i >= 0; i--) {
        for (int j = 0; j < energyMap[i].size(); j++) {
            if (j == 0) {
                // process leftmost index
                cumulativeEnergyMatrix[i][j] = energyMap[i][j] + min(cumulativeEnergyMatrix[i+1][j], cumulativeEnergyMatrix[i+1][j+1]);
            } else if (j == energyMap[i].size()-1) {
                // process rightmost index
                cumulativeEnergyMatrix[i][j] = energyMap[i][j] + min(cumulativeEnergyMatrix[i+1][j], cumulativeEnergyMatrix[i+1][j-1]);
            } else {
                // process everything else
                cumulativeEnergyMatrix[i][j] = energyMap[i][j] + min(cumulativeEnergyMatrix[i+1][j], min(cumulativeEnergyMatrix[i+1][j-1], cumulativeEnergyMatrix[i+1][j+1]));
            }
        }
    }

    return cumulativeEnergyMatrix;
}

// REMOVE SEAM FUNCTION
// input: cumulative energy matrix
// output: vector holding the path to traverse for seam removal
vector<int> findSeam(vector<vector<int>> cumulativeEnergyMatrix) {

    vector<int> seam(cumulativeEnergyMatrix.size());
    // get lowest value in first row of cumulative energy matrix
    int lowestValue = cumulativeEnergyMatrix[0][0], lowestValueIndex = 0; 
    for (int j = 0; j < cumulativeEnergyMatrix[0].size(); j++) {
        if (cumulativeEnergyMatrix[0][j] < lowestValue) {
            lowestValue = cumulativeEnergyMatrix[0][j];
            lowestValueIndex = j;
        }
    }

    // add this value to our seam vector
    seam[0] = lowestValueIndex;

    // by this point we should know where to begin carving from. This is our lowestValueIndex variable.
    for (int i = 1; i < cumulativeEnergyMatrix.size(); i++) {
        
        if (seam[i-1] == 0) {
            lowestValue = min(cumulativeEnergyMatrix[i][seam[i-1]], cumulativeEnergyMatrix[i][seam[i-1]+1]);
        } else if (seam[i-1] == cumulativeEnergyMatrix[i].size()-1) {
            lowestValue = min(cumulativeEnergyMatrix[i][seam[i-1]], cumulativeEnergyMatrix[i][seam[i-1]-1]);
        } else {
            lowestValue = min(cumulativeEnergyMatrix[i][seam[i-1]], min(cumulativeEnergyMatrix[i][seam[i-1]-1], cumulativeEnergyMatrix[i][seam[i-1]+1]));
        }

        if (seam[i-1] != cumulativeEnergyMatrix[i].size()-1) {
            if (lowestValue == cumulativeEnergyMatrix[i][seam[i-1]+1]) {
                lowestValueIndex = seam[i-1]+1;
            } 
        }

        if (lowestValue == cumulativeEnergyMatrix[i][seam[i-1]]) {
            lowestValueIndex = seam[i-1];
        } 

        if (seam[i-1] != 0) {
            if (lowestValue == cumulativeEnergyMatrix[i][seam[i-1]-1]) {
                lowestValueIndex = seam[i-1]-1;
            }
        }

        seam[i] = lowestValueIndex;
    }

    return seam;
}

// REMOVE SEAMS FUNCTION
// input: value left of index, value right of index, value above index, value below index, value of index, selector
// output: energy at index
void removeSeams(vector<vector<int>> &matrix, vector<int> seam) {
    
    for (int i = 0; i < matrix.size(); i++) {
        for (auto it = matrix[i].begin(); it < matrix[i].end(); it++) {
            int index = it - matrix[i].begin();
            if (index == seam[i]) {
                it = matrix[i].erase(it);
            }
        }
    }
}

// REMOVE SEAMS DRIVER FUNCTION
// input: matrix, number of vertical columns to remove, number of horizontal rows to remove
// output: matrix
vector<vector<int>> removeSeamsDriver(vector<vector<int>> mat, int numV, int numH) {

    vector<vector<int>> originalMatrix(mat.size(), vector<int>(mat[0].size()));
    originalMatrix = mat;

    vector<vector<int>> energyMatrix(originalMatrix.size(), vector<int>(originalMatrix[0].size()));
    vector<vector<int>> cumulativeEnergyMatrix(energyMatrix.size(), vector<int>(energyMatrix[0].size()));

    if (numV > 0) {
        for (int v = 0; v < numV; v++) {
            energyMatrix = createEnergyMatrix(originalMatrix);

            cumulativeEnergyMatrix = createCumulativeEnergyMatrix(energyMatrix);

            vector<int> seam(cumulativeEnergyMatrix.size());
            seam = findSeam(cumulativeEnergyMatrix);

            removeSeams(originalMatrix, seam);
            removeSeams(energyMatrix, seam);
            removeSeams(cumulativeEnergyMatrix, seam);
        }
    }

    vector<vector<int>> transposedMatrix(originalMatrix[0].size(), vector<int>(originalMatrix.size()));
    if (numH > 0) transposedMatrix = transposeMatrix(originalMatrix);

    vector<vector<int>> transposedEnergyMatrix(transposedMatrix.size(), vector<int>(transposedMatrix[0].size()));
    vector<vector<int>> transposedCumulativeEnergyMatrix(transposedEnergyMatrix.size(), vector<int>(transposedEnergyMatrix[0].size()));

    if (numH > 0) {
        for (int h = 0; h < numH; h++) {
            transposedEnergyMatrix = createEnergyMatrix(transposedMatrix);

            transposedCumulativeEnergyMatrix = createCumulativeEnergyMatrix(transposedEnergyMatrix);

            vector<int> seam(transposedCumulativeEnergyMatrix.size());
            seam = findSeam(transposedCumulativeEnergyMatrix);

            removeSeams(transposedMatrix, seam);
            removeSeams(transposedEnergyMatrix, seam);
            removeSeams(transposedCumulativeEnergyMatrix, seam);
        }
    }
    
    vector<vector<int>> finalMatrix(transposedMatrix[0].size(), vector<int>(transposedMatrix.size()));
    if (numH > 0) finalMatrix = transposeMatrix(transposedMatrix);

    if (numH > 0) {
        return finalMatrix;
    }

    return originalMatrix;
}