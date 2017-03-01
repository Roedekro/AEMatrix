#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>

using namespace std;

int* buildSquareMatrix(int m, int range) {

    int* matrix = new int[m*m];
    for(int i = 0; i < m*m; i++) {
        matrix[i] = (rand() % range) + 1;
    }
    return matrix;
}

int* buildSkewedMatrix(int m, int n, int range) {

    int* matrix = new int[m*n];
    for(int i = 0; i < m*n; i++) {
        matrix[i] = (rand() % range) + 1;
    }
    return matrix;
}

int* multiplyMatrixRow(int* a, int m, int n, int* b, int p) {

    int* c = new int[m*p];
    for(int i = 0; i < m*p; i++) {
        c[i] = 0;
    }

    // matrix[i][j]
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < p; j++) {
            for(int k = 0; k < n; k++) {
                // c[i][j] += a[i][k] * b[k][j]
                c[i*n+j] += a[i*n+k] * b[k*n+j];
            }
        }
    }

    return c;

}

int* multiplyMatrixColumn(int* a, int m, int n, int* b, int p) {

    int* c = new int[m*p];
    for(int i = 0; i < m*p; i++) {
        c[i] = 0;
    }

    // matrix[i][j]
    for(int j = 0; j < p; j++) {
        for(int i = 0; i < p; i++) {
            for(int k = 0; k < n; k++) {
                // c[i][j] += a[i][k] * b[k][j]
                c[j*n+i] += a[k*n+i] * b[j*n+k];
            }
        }
    }

    return c;

}

/*
 * A is row ordered, B is column ordered
 */

int* multiplyMatrixMixed(int* a, int m, int n, int* b, int p) {

    int* c = new int[m*p];
    for(int i = 0; i < m*p; i++) {
        c[i] = 0;
    }

    // matrix[i][j]
    for(int i = 0; i < m; i++) {
        for(int j = 0; j < p; j++) {
            for(int k = 0; k < n; k++) {
                // c[i][j] += a[i][k] * b[k][j]
                c[i*n+j] += a[i*n+k] * b[j*n+k];
            }
        }
    }

    return c;

}

//void testRowVsColumn(int m, int n, int p, int range) {

void testRowVsColumnVsMixedSquare(int s, int offset, int increment, int range) {

    int* timeRow = new int[(s-offset)/increment];
    int* timeCol = new int[(s-offset)/increment];
    int* timeMix = new int[(s-offset)/increment];
    int counter = -1;

    for(int x = offset; x <= s; x=x+increment) {

        counter++;
        int m = x;
        int n = x;
        int p = x;

        cerr << "Building " << x << '\n';
        int* a = buildSkewedMatrix(m,n,range);
        int* b = buildSkewedMatrix(n,p,range);
        cerr << "Finished\n";

        typedef std::chrono::system_clock Clock;
        auto start = Clock::now();

        int* cRow = multiplyMatrixRow(a,m,n,b,p);

        auto stop = Clock::now();
        auto total = stop-start;
        long millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeRow[counter] = millis;
        cerr << cRow[1] << '\n';

        int* cRow2 = multiplyMatrixRow(a,m,n,b,p);

        cerr << cRow2[1] << '\n';


        start = Clock::now();

        int* cCol = multiplyMatrixColumn(a,m,n,b,p);

        stop = Clock::now();
        total = stop-start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeCol[counter] = millis;
        cerr << cCol[1] << '\n';

        int* cCol2 = multiplyMatrixColumn(a,m,n,b,p);

        cerr << cCol2[1] << '\n';


        start = Clock::now();

        int* cMix = multiplyMatrixMixed(a,m,n,b,p);

        stop = Clock::now();
        total = stop-start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeMix[counter] = millis;
        cerr << cMix[1] << '\n';

        int* cMix2 = multiplyMatrixMixed(a,m,n,b,p);

        cerr << cMix2[1] << '\n';
    }

    for(int i = 0; i <= counter; i++) {
        cout << (offset+i*increment) << '\t' << timeRow[i] << '\t' << timeCol[i] << '\t' << timeMix[i] << '\n';
    }



}



int main(int argc, char* argv[]) {

    int size,test,range,offset,increment;
    if(argc != 5) {
        cout << "Arguments are <test> <max size> <min size> <increment> <range>\n";
        test = 1;
        size = 700;
        offset = 500;
        increment = 100;
        range = 100000;
    }
    else {
        test = atoi(argv[1]);
        size = atoi(argv[2]);
        offset = atoi(argv[3]);
        increment = atoi(argv[4]);
        range = atoi(argv[5]);
    }


    if(test == 1) {
        testRowVsColumnVsMixedSquare(size,offset,increment,range);
    }

    /*int* test1 = new int[4]{1,2,3,4};
    int* test2 = new int[4]{1,2,3,4};

    int* res1 = multiplyMatrixRow(test1,2,2,test2,2);
    int* res2 = multiplyMatrixColumn(test1,2,2,test2,2);

    for(int i = 0; i < 4; i++) {
        cout << res1[i] << '\n';
    }

    cout << "---\n";

    for(int i = 0; i < 4; i++) {
        cout << res2[i] << '\n';
    }*/

    return 0;
}