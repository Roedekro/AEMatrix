#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>

using namespace std;

int* transposeNaive(int* a, int* b, int m, int n) {

    for(int i = 0; i < m; i++) {
        for(int j = 0; j < n; j++) {
            b[j*m+i] = a[i*n+j];
        }
    }
}

void transposeRec(int* a, int* b, int mStart, int nStart, int mStop, int nStop, int origM, int origN) {

    int newN, newM;
    int mDif = mStop - mStart;
    int nDif = nStop - nStart;
    if(mDif == 1 || nDif == 1) {
        /*cout << "Transposing " << mStart << "," << nStart << "," << mStop << "," << nStop
                << "," << origM << "," << origN << '\n';*/
        for(int i = 0; i < mDif; i++) {
            for(int j = 0; j < nDif; j++) {
                b[(j+nStart)*origM+mStart+i] = a[mStart*origN+i*origN+nStart+j];
                /*cout << "Placed " << a[mStart*origN+i*origN+nStart+j] << " from " << mStart*origN+i*origN+nStart+j
                    << " to " << (j+nStart)*origM+mStart+i << '\n';*/
            }
        }
    }
    else {

        if(nDif >= mDif) {
            // Split langs n
            int lengthN = nDif/2;
            newN = nStart+lengthN;
            transposeRec(a,b,mStart,nStart,mStop,newN,origM,origN);
            transposeRec(a,b,mStart,newN,mStop,nStop,origM,origN);
        }
        else {
            // Split langs m
            int lengthM = mDif/2;
            newM = mStart + lengthM;
            transposeRec(a,b,mStart,nStart,newM,nStop,origM,origN);
            transposeRec(a,b,newM,nStart,mStop,nStop,origM,origN);
        }
    }
}

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

        delete[] a;
        delete[] b;
        delete[] cRow;
        delete[] cRow2;
        delete[] cCol;
        delete[] cCol2;
        delete[] cMix;
        delete[] cMix2;
    }

    for(int i = 0; i <= counter; i++) {
        cout << (offset+i*increment) << '\t' << timeRow[i] << '\t' << timeCol[i] << '\t' << timeMix[i] << '\n';
    }



}

void testTranspose(int s, int offset, int increment, int range) {

    int* timeClassic = new int[(s-offset)/increment];
    int* timeRecursive = new int[(s-offset)/increment];
    int counter = -1;

    for(int x = offset; x <= s; x=x+increment) {
        counter++;
        int m = x;
        int n = x;
        int p = x;

        cerr << "Building " << x << '\n';
        int* a = buildSkewedMatrix(m,n,range);
        cerr << "Finished\n";
        int* b = new int[m*n];

        typedef std::chrono::system_clock Clock;
        auto start = Clock::now();

        transposeNaive(a,b,m,n);

        auto stop = Clock::now();
        auto total = stop-start;
        long millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeClassic[counter] = millis;
        cerr << b[1] << '\n';

        transposeNaive(a,b,m,n);

        cerr << b[1] << '\n';


        start = Clock::now();

        transposeRec(a,b,0,0,m,n,m,n);

        stop = Clock::now();
        total = stop-start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeRecursive[counter] = millis;
        cerr << b[1] << '\n';

        transposeRec(a,b,0,0,m,n,m,n);

        cerr << b[1] << '\n';



        delete[] a;
        delete[] b;
    }

    for(int i = 0; i <= counter; i++) {
        cout << (offset+i*increment) << '\t' << timeClassic[i] << '\t' << timeRecursive[i] << '\n';
    }

}

void testMultiplyTranspose (int s, int offset, int increment, int range) {

    int* timeNoTranspose = new int[(s-offset)/increment];
    int* timeWithTranspose = new int[(s-offset)/increment];
    int counter = -1;

    for(int x = offset; x <= s; x=x+increment) {
        counter++;
        int m = x;
        int n = x;
        int p = x;

        cerr << "Building " << x << '\n';
        int* a = buildSkewedMatrix(m,n,range);
        int* b = buildSkewedMatrix(m,n,range);
        cerr << "Finished\n";
        int* c = new int[m*n];

        typedef std::chrono::system_clock Clock;
        auto start = Clock::now();

        int* cRow = multiplyMatrixRow(a,m,n,c,p);

        auto stop = Clock::now();
        auto total = stop-start;
        long millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeNoTranspose[counter] = millis;
        cerr << cRow[1] << '\n';

        int* cRow2 = multiplyMatrixRow(a,m,n,c,p);

        cerr << cRow2[1] << '\n';


        start = Clock::now();

        transposeRec(b,c,0,0,m,n,m,n);
        int* cMix = multiplyMatrixMixed(a,m,n,c,p);

        stop = Clock::now();
        total = stop-start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeWithTranspose[counter] = millis;
        cerr << cMix[1] << '\n';

        transposeRec(b,c,0,0,m,n,m,n);
        int* cMix2 = multiplyMatrixMixed(a,m,n,c,p);

        cerr << cMix2[1] << '\n';



        delete[] a;
        delete[] b;
        delete[] c;
        delete[] cRow;
        delete[] cRow2;
        delete[] cMix;
        delete[] cMix2;
    }

    for(int i = 0; i <= counter; i++) {
        cout << (offset+i*increment) << '\t' << timeNoTranspose[i] << '\t' << timeWithTranspose[i] << '\n';
    }

}



int main(int argc, char* argv[]) {

    int size,test,range,offset,increment;
    if(argc != 5) {
        cout << "Arguments are <test> <max size> <min size> <increment> <range>\n";
        test = 3;
        size = 1000;
        offset = 1000;
        increment = 1000;
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
    else if(test == 2) {
        testTranspose(size,offset,increment,range);
    }
    else if(test == 3) {
        testMultiplyTranspose(size,offset,increment,range);
    }

    // Test transpose functions
    /*int* test1 = new int[9]{1,2,3,4,5,6,7,8,9};
    int* test2 = new int[9];
    int* test3 = new int[9];

    transposeNaive(test1,test2,3,3);
    transposeRec(test1,test3,0,0,3,3,3,3);

    for(int i = 0; i < 9; i++) {
        cout << test2[i] << '\n';
    }

    cout << "---\n";

    for(int i = 0; i < 9; i++) {
        cout << test3[i] << '\n';
    }*/

    // Test multiply functions
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