#include <iostream>
#include <stdlib.h>
#include <ctime>
#include <chrono>
#include <math.h>
#include <algorithm>

using namespace std;

void plus(int* a, int* b, int* c, int n) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            c[i*n+j] = a[i*n+j] + b[i*n+j];
        }
    }
}

void minus(int* a, int* b, int* c, int n) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            c[i*n+j] = a[i*n+j] + b[i*n+j];
        }
    }
}

// Virker kun med n=2^x
// Se https://en.wikipedia.org/wiki/Strassen_algorithm
void strassen(int* a, int* b, int* c, int n) {

    if(n == 1) {
        c[0] = a[0] * b[0];
        return;
    }

    int newN = n/2;

    int* a11 = new int[newN*newN];
    int* a12 = new int[newN*newN];
    int* a21 = new int[newN*newN];
    int* a22 = new int[newN*newN];

    int* b11 = new int[newN*newN];
    int* b12 = new int[newN*newN];
    int* b21 = new int[newN*newN];
    int* b22 = new int[newN*newN];

    int* c11 = new int[newN*newN];
    int* c12 = new int[newN*newN];
    int* c21 = new int[newN*newN];
    int* c22 = new int[newN*newN];

    int* m1 = new int[newN*newN];
    int* m2 = new int[newN*newN];
    int* m3 = new int[newN*newN];
    int* m4 = new int[newN*newN];
    int* m5 = new int[newN*newN];
    int* m6 = new int[newN*newN];
    int* m7 = new int[newN*newN];

    int* resA = new int[newN*newN];
    int* resB = new int[newN*newN];

    // Fyld a11 .. a22 og b11 .. b22 ud
    for(int i = 0; i < newN; i++) {
        for(int j = 0; j < newN; j++) {
            a11[i*newN+j] = a[i*n+j];
        }
    }
    for(int i = 0; i < newN; i++) {
        for(int j = newN; j < n; j++) {
            a12[i*newN+j] = a[i*n+j];
        }
    }

    for(int i = newN; i < n; i++) {
        for(int j = 0; j < newN; j++) {
            a21[i*newN+j] = a[i*n+j];
        }
    }

    for(int i = newN; i < n; i++) {
        for(int j = newN; j < n; j++) {
            a22[i*newN+j] = a[i*n+j];
        }
    }

    // B

    for(int i = 0; i < newN; i++) {
        for(int j = 0; j < newN; j++) {
            b11[i*newN+j] = b[i*n+j];
        }
    }
    for(int i = 0; i < newN; i++) {
        for(int j = newN; j < n; j++) {
            b12[i*newN+j] = b[i*n+j];
        }
    }

    for(int i = newN; i < n; i++) {
        for(int j = 0; j < newN; j++) {
            b21[i*newN+j] = b[i*n+j];
        }
    }

    for(int i = newN; i < n; i++) {
        for(int j = newN; j < n; j++) {
            b22[i*newN+j] = b[i*n+j];
        }
    }

    // Udregn M matricerne

    // M1 = (a11 + a22) * (b11 + b22)
    plus(a11,a22,resA,newN);
    plus(b11,b22,resB,newN);
    strassen(resA,resB,m1,newN);

    // M2 = (a21 + a22) * b11
    plus(a21,a22,resA,newN);
    strassen(resA,b11,m2,newN);

    // M3 = a11 * (b12 - b22)
    minus(b12,b22,resB,newN);
    strassen(a11,resB,m3,newN);

    // M4 = a22 * (b21 - b11)
    minus(b21,b11,resB,newN);
    strassen(a22,resB,m4,newN);

    // M5 = (a11 + a12) * b22
    plus(a11,a12,resA,newN);
    strassen(resA,b22,m5,newN);

    // M6 = (a21 - a11) * (b11 + b12)
    minus(a21,a11,resA,newN);
    plus(b11,b12,resB,newN);
    strassen(resA,resB,m6,newN);

    // M7 = (a12 - a22) * (b21 + b22)
    minus(a12,a22,resA,newN);
    plus(b21,b22,resB,newN);
    strassen(resA,resB,m7,newN);

    // Så finder vi C

    // c11 = m1 + m4 - m5 + m7
    plus(m1,m4,resA,newN);
    minus(resA,m5,resB,newN);
    plus(resB,m7,c11,newN);

    // c12 = m3 + m5
    plus(m3,m5,c12,newN);

    // c21 = m2 + m4
    plus(m2,m4,c21,newN);

    // c22 = m1 - m2 + m3 + m6
    minus(m1,m2,resA,newN);
    plus(resA,m3,resB,newN);
    plus(resB,m6,c22,newN);

    // Fyld c11 ind i C
    for(int i = 0; i < newN; i++) {
        for(int j = 0; j < newN; j++) {
            c[i*n+j] = c11[i*newN+j];
        }
    }

    // Fyld c12 ind i C
    for(int i = 0; i < newN; i++) {
        for(int j = newN; j < n; j++) {
            c[i*n+j] = c11[i*newN+j];
        }
    }

    // Fyld c21 ind i C
    for(int i = newN; i < n; i++) {
        for(int j = 0; j < newN; j++) {
            c[i*n+j] = c11[i*newN+j];
        }
    }

    // Fyld c22 ind i C
    for(int i = newN; i < n; i++) {
        for(int j = newN; j < n; j++) {
            c[i*n+j] = c11[i*newN+j];
        }
    }

    // Ryd pænt op

    delete[] a11;
    delete[] a12;
    delete[] a21;
    delete[] a22;

    delete[] b11;
    delete[] b12;
    delete[] b21;
    delete[] b22;

    delete[] c11;
    delete[] c12;
    delete[] c21;
    delete[] c22;

    delete[] m1;
    delete[] m2;
    delete[] m3;
    delete[] m4;
    delete[] m5;
    delete[] m6;
    delete[] m7;

    delete[] resA;
    delete[] resB;

}

// Virker kun med n = 2^x
int* multiplyRecursiveLayoutSquare(int* a, int* b, int* c, int n) {

    // matrix[i][j]
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            for(int k = 0; k < n; k++) {

                // Figure out positions
                unsigned short i2 = (short) i;
                unsigned short j2 = (short) j;
                unsigned short k2 = (short) k;
                unsigned int pos1 = 0;
                unsigned int pos2 = 0;
                unsigned int pos3 = 0;
                // z |= (x & 1U << i) << i | (y & 1U << i) << (i + 1);
                // Taken from http://graphics.stanford.edu/~seander/bithacks.html#InterleaveTableObvious
                for (int x = 0; x < 16; x++) {
                    pos1 |= (i2 & 1U << x) << x | (j2 & 1U << x) << (x + 1);
                }
                for (int x = 0; x < 16; x++) {
                    pos2 |= (i2 & 1U << x) << x | (k2 & 1U << x) << (x + 1);
                }
                for (int x = 0; x < 16; x++) {
                    pos3 |= (k2 & 1U << x) << x | (j2 & 1U << x) << (x + 1);
                }

                // c[i][j] += a[i][k] * b[k][j]
                c[pos1] += a[pos2] * b[pos3];

                cout << i << "," << j << "," << k << '\n';
                cout << pos1 << "," << pos2 << "," << pos3 << '\n';
                cout << a[pos2] * b[pos3] << "," << a[pos2] << "," << b[pos3] << '\n';
                cout << "---\n";
            }
        }
    }

    return c;

}

int* multiplyTiledSquare(int* a, int* b, int* c, int n, int s) {

    int to = n/s;
    if(n % s != 0) {
        to++;
    }

    for(int x = 0; x < to; x++) {
        for(int y = 0; y < to; y++) {
            for(int z = 0; z < to; z++) {
                // matrix[i][j]
                for(int i = x*s; (i < (x+1)*s && i < n); i++) {
                    for(int j = y*s; (j < (y+1)*s && j < n); j++) {
                        for(int k = z*s; (k < (z+1)*s && k < n); k++) {
                            // c[i][j] += a[i][k] * b[k][j]
                            c[i*n+j] += a[i*n+k] * b[k*n+j];
                            cout << x << "," << y << "," << z << "\t" << i << "," << j << "," << z << '\n';
                        }
                    }
                }
            }
        }

    }
    return c;
}

void multiplyRecursiveMixed(int* a, int* b, int* c, int mStart, int nStart, int pStart, int mStop, int nStop, int pStop,
                             int mOrig, int nOrig, int pOrig) {

    int mDif = mStop-mStart;
    int nDif = nStop-nStart;
    int pDif = pStop-pStart;
    //cout << mStart << "," << nStart << "," << pStart << "\t" << mDif << "," << nDif << "," << pDif << '\n';
    //if(mDif == 1 && pDif == 1) {
    /*if(mDif == 1 || pDif == 1) {
        for(int i = mStart; i < mStop; i++) {
            for(int j = pStart; j < pStop; j++) {
                for(int k = nStart; k < nStop; k++) {
                    // Mixed c[i*n+j] += a[i*n+k] * b[j*n+k];
                    c[i*pOrig+j] += a[i*nOrig+k] * b[j*pOrig+k];
                }
            }
        }
    }*/
    if(mDif == 1 && nDif == 1 && pDif == 1) {
        c[mStart*pOrig+pStart] += a[mStart*nOrig+nStart] * b[pStart*pOrig+nStart];
        /*cout << "Added " << a[mStart*nOrig+nStart] * b[pStart*pOrig+nStart] << " to " << mStart*pOrig+pStart << " totalling "
        << c[mStart*pOrig+pStart] << '\n';*/
    }
    else {
        if(mDif >= max(nDif, pDif)) {
            int newM = (mStart + mStop)/2;
            multiplyRecursiveMixed(a,b,c,mStart,nStart,pStart,newM,nStop,pStop,mOrig,nOrig,pOrig);
            multiplyRecursiveMixed(a,b,c,newM,nStart,pStart,mStop,nStop,pStop,mOrig,nOrig,pOrig);
        }
        else if(nDif >= max(mDif,pDif)) {
            int newN = (nStart+nStop)/2;
            multiplyRecursiveMixed(a,b,c,mStart,nStart,pStart,mStop,newN,pStop,mOrig,nOrig,pOrig);
            multiplyRecursiveMixed(a,b,c,mStart,newN,pStart,mStop,nStop,pStop,mOrig,nOrig,pOrig);
        }
        else {
            int newP = (pStart+pStop)/2;
            multiplyRecursiveMixed(a,b,c,mStart,nStart,pStart,mStop,nStop,newP,mOrig,nOrig,pOrig);
            multiplyRecursiveMixed(a,b,c,mStart,nStart,newP,mStop,nStop,pStop,mOrig,nOrig,pOrig);
        }
    }
}

void multiplyRecursiveRow(int* a, int* b, int* c, int mStart, int nStart, int pStart, int mStop, int nStop, int pStop,
                            int mOrig, int nOrig, int pOrig) {

    int mDif = mStop-mStart;
    int nDif = nStop-nStart;
    int pDif = pStop-pStart;
    /*if(mDif == 1 || pDif == 1) {
        for(int i = mStart; i < mStop; i++) {
            for(int j = pStart; j < pStop; j++) {
                for(int k = nStart; k < nStop; k++) {
                    // Row c[i*n+j] += a[i*n+k] * b[k*n+j];
                    c[i*pOrig+j] += a[i*nOrig+k] * b[k*pOrig+j];
                }
            }
        }
    }*/
    //cout << mStart << "," << nStart << "," << pStart << "\t" << mDif << "," << nDif << "," << pDif << '\n';
    if(mDif == 1 && nDif == 1 && pDif == 1) {
        c[mStart*pOrig+pStart] += a[mStart*nOrig+nStart] * b[nStart*pOrig+pStart];
        /*cout << "Added " << a[mStart*nOrig+nStart] * b[nStart*pOrig+pStart] << " to " << mStart*pOrig+pStart << " totalling "
             << c[mStart*pOrig+pStart] << '\n';*/
    }
    else {
        if(mDif >= max(nDif, pDif)) {
            int newM = (mStart + mStop)/2;
            multiplyRecursiveRow(a,b,c,mStart,nStart,pStart,newM,nStop,pStop,mOrig,nOrig,pOrig);
            multiplyRecursiveRow(a,b,c,newM,nStart,pStart,mStop,nStop,pStop,mOrig,nOrig,pOrig);
        }
        else if(nDif >= max(mDif,pDif)) {
            int newN = (nStart+nStop)/2;
            multiplyRecursiveRow(a,b,c,mStart,nStart,pStart,mStop,newN,pStop,mOrig,nOrig,pOrig);
            multiplyRecursiveRow(a,b,c,mStart,newN,pStart,mStop,nStop,pStop,mOrig,nOrig,pOrig);
        }
        else {
            int newP = (pStart+pStop)/2;
            multiplyRecursiveRow(a,b,c,mStart,nStart,pStart,mStop,nStop,newP,mOrig,nOrig,pOrig);
            multiplyRecursiveRow(a,b,c,mStart,nStart,newP,mStop,nStop,pStop,mOrig,nOrig,pOrig);
        }
    }
}


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
    //if(mDif == 1 || nDif == 1) {
    if(mDif == 1 && nDif == 1) {
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

int* multiplyMatrixRow(int* a, int m, int n, int* b, int p, int* c) {

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

int* multiplyMatrixColumn(int* a, int m, int n, int* b, int p, int* c) {

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

int* multiplyMatrixMixed(int* a, int m, int n, int* b, int p, int* c) {

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

//void testRowVsColumnVsMixedSquare(int s, int offset, int increment, int range) {

void testRowVsColumnVsMixedSquare(int s, int runs, int range) {

    /*int* timeRow = new int[(s-offset)/increment+1];
    int* timeCol = new int[(s-offset)/increment+1];
    int* timeMix = new int[(s-offset)/increment+1];*/
    int* timeRow = new int[s-1];
    int* timeCol = new int[s-1];
    int* timeMix = new int[s-1];

    int counter = -1;

    for(int x = 2; x <= s; x++) {

        int temp = pow(3,x);
        counter++;
        int m = temp;
        int n = temp;
        int p = temp;

        cerr << "Building " << x << '\n';
        int* a = buildSkewedMatrix(m,n,range);
        int* b = buildSkewedMatrix(n,p,range);
        int* c = new int[m*p];
        cerr << "Finished\n";

        typedef std::chrono::system_clock Clock;
        auto start = Clock::now();

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixRow(a,m,n,b,p,c);
        }



        auto stop = Clock::now();
        auto total = stop-start;
        long millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeRow[counter] = millis;
        cerr << c[1] << '\n';

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixRow(a,m,n,b,p,c);
        }

        cerr << c[1] << '\n';


        start = Clock::now();

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixColumn(a,m,n,b,p,c);
        }

        stop = Clock::now();
        total = stop-start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeCol[counter] = millis;
        cerr << c[1] << '\n';

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixColumn(a,m,n,b,p,c);
        }

        cerr << c[1] << '\n';


        start = Clock::now();

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixMixed(a,m,n,b,p,c);
        }

        stop = Clock::now();
        total = stop-start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeMix[counter] = millis;
        cerr << c[1] << '\n';

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixMixed(a,m,n,b,p,c);
        }

        cerr << c[1] << '\n';

        delete[] a;
        delete[] b;
        delete[] c;
    }

    for(int i = 0; i <= counter; i++) {
        cout << (i+2) << '\t' << timeRow[i] << '\t' << timeCol[i] << '\t' << timeMix[i] << '\n';
    }



}

void testTranspose(int s, int runs, int range) {

    int* timeClassic = new int[s-1];
    int* timeRecursive = new int[s-1];
    int counter = -1;

    for(int x = 2; x <= s; x++) {
        int temp = pow(3,x);
        counter++;
        int m = temp;
        int n = temp;
        int p = temp;

        cerr << "Building " << x << '\n';
        int* a = buildSkewedMatrix(m,n,range);
        int* b = new int[m*n];
        cerr << "Finished\n";


        typedef std::chrono::system_clock Clock;
        auto start = Clock::now();

        for(int j = 0; j < runs; j++) {
            transposeNaive(a,b,m,n);
        }


        auto stop = Clock::now();
        auto total = stop-start;
        long millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeClassic[counter] = millis;
        cerr << b[1] << '\n';

        for(int j = 0; j < runs; j++) {
            transposeNaive(a,b,m,n);
        }

        cerr << b[1] << '\n';


        start = Clock::now();

        for(int j = 0; j < runs; j++) {
            transposeRec(a,b,0,0,m,n,m,n);
        }


        stop = Clock::now();
        total = stop-start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeRecursive[counter] = millis;
        cerr << b[1] << '\n';

        for(int j = 0; j < runs; j++) {
            transposeRec(a,b,0,0,m,n,m,n);
        }

        cerr << b[1] << '\n';



        delete[] a;
        delete[] b;
    }

    for(int i = 0; i <= counter; i++) {
        cout << (i+2) << '\t' << timeClassic[i] << '\t' << timeRecursive[i] << '\n';
    }

}

void testMultiplyTranspose (int s, int runs, int range) {

    int* timeNoTranspose = new int[s-1];
    int* timeWithTranspose = new int[s-1];
    int counter = -1;

    for(int x = 2; x <= s; x++) {
        int temp = pow(3,x);
        counter++;
        int m = temp;
        int n = temp;
        int p = temp;

        cerr << "Building " << x << '\n';
        int* a = buildSkewedMatrix(m,n,range);
        int* b = buildSkewedMatrix(m,n,range);
        int* c = new int[m*p];
        int* d = new int[m*n];
        cerr << "Finished\n";


        typedef std::chrono::system_clock Clock;
        auto start = Clock::now();

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixRow(a,m,n,b,p,c);
        }

        auto stop = Clock::now();
        auto total = stop-start;
        long millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeNoTranspose[counter] = millis;
        cerr << c[1] << '\n';

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixRow(a,m,n,b,p,c);
        }

        cerr << c[1] << '\n';


        start = Clock::now();

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            transposeNaive(b,d,m,n);
            multiplyMatrixMixed(a,m,n,d,p,c);
        }

        stop = Clock::now();
        total = stop-start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeWithTranspose[counter] = millis;
        cerr << c[1] << '\n';

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            transposeNaive(b,d,m,n);
            multiplyMatrixMixed(a,m,n,d,p,c);
        }

        cerr << c[1] << '\n';



        delete[] a;
        delete[] b;
        delete[] c;
        delete[] d;
    }

    for(int i = 0; i <= counter; i++) {
        cout << (i+2) << '\t' << timeNoTranspose[i] << '\t' << timeWithTranspose[i] << '\n';
    }

}

void testMultiplyRecursive (int s, int runs, int range, int realS) {

    int *timeNormal = new int[s-1];
    int *timeRecursive = new int[s-1];
    int *timeTiled = new int[s-1];
    int counter = -1;

    for(int x = 2; x <= s; x++) {
        int temp = pow(3,x);
        counter++;
        int m = temp;
        int n = temp;
        int p = temp;

        cerr << "Building " << x << '\n';
        int *a = buildSkewedMatrix(m, n, range);
        int *b = buildSkewedMatrix(n, p, range);
        int *c = new int[m * p];
        cerr << "Finished\n";



        typedef std::chrono::system_clock Clock;
        auto start = Clock::now();

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixRow(a,m,n,b,p,c);
        }

        auto stop = Clock::now();
        auto total = stop - start;
        long millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeNormal[counter] = millis;
        cerr << c[1] << '\n';

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyMatrixRow(a,m,n,b,p,c);
        }

        cerr << c[1] << '\n';


        start = Clock::now();

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyRecursiveRow(a, b, c, 0, 0, 0, m, n, p, m, n, p);
        }

        stop = Clock::now();
        total = stop - start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeRecursive[counter] = millis;
        cerr << c[1] << '\n';

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyRecursiveRow(a, b, c, 0, 0, 0, m, n, p, m, n, p);
        }

        cerr << c[1] << '\n';


        start = Clock::now();

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyTiledSquare(a,b,c,n,realS);
        }

        stop = Clock::now();
        total = stop - start;
        millis = std::chrono::duration_cast<std::chrono::milliseconds>(total).count();
        timeTiled[counter] = millis;
        cerr << c[1] << '\n';

        for(int j = 0; j < runs; j++) {
            for(int i = 0; i < m*p; i++) {
                c[i] = 0;
            }
            multiplyTiledSquare(a,b,c,n,realS);
        }

        cerr << c[1] << '\n';


        delete[] a;
        delete[] b;
        delete[] c;

    }
    for (int i = 0; i <= counter; i++) {
        cout << (i+2) << '\t' << timeNormal[i] << '\t' << timeRecursive[i] << '\t' << timeTiled[i] << '\n';
    }
}



int main(int argc, char* argv[]) {

    int size,test,range,offset,increment,runs,s;
    if(argc != 6) {
        cout << "Arguments are <test> <pow> <range> <runs> <s>\n";
        test = 1;
        size = 6;
        offset = 10000;
        increment = 1000;
        range = 100000;
        runs = 1000;
        s = 4;
    }
    else {
        test = atoi(argv[1]);
        size = atoi(argv[2]);
        //offset = atoi(argv[3]);
        //increment = atoi(argv[4]);
        range = atoi(argv[3]);
        runs = atoi(argv[4]);
        s = atoi(argv[5]);

    }


    if(test == 1) {
        testRowVsColumnVsMixedSquare(size,runs,range);
    }
    else if(test == 2) {
        testTranspose(size,runs,range);
    }
    else if(test == 3) {
        testMultiplyTranspose(size,runs,range);
    }
    else if(test == 4) {
        testMultiplyRecursive(size,runs,range,s);
    }

    /* Correct result will be
    30
    36
    42
    66
    81
    96
    102
    126
    150
     */


    // Test multiplyRecursiveStructureSquare
    /*int* test1 = new int[16]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    int* test2 = new int[16]{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
    int* test3 = new int[16]{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    multiplyRecursiveLayoutSquare(test1,test2,test3,3);

    for(int i = 0; i < 9; i++) {
        cout << test3[i] << '\n';
    }*/

    // Test tiledMultSquare
    /*int* test1 = new int[9]{1,2,3,4,5,6,7,8,9};
    int* test2 = new int[9]{1,2,3,4,5,6,7,8,9};
    int* test3 = new int[9]{0,0,0,0,0,0,0,0,0};

    multiplyTiledSquare(test1,test2,test3,3,3);

    for(int i = 0; i < 9; i++) {
        cout << test3[i] << '\n';
    }*/


    // Test recursive multiply
    /*int* test1 = new int[9]{1,2,3,4,5,6,7,8,9};
    int* test2 = new int[9]{1,2,3,4,5,6,7,8,9};
    int* test3 = new int[9]{0,0,0,0,0,0,0,0,0};

    multiplyRecursiveRow(test1,test2,test3,0,0,0,3,3,3,3,3,3);

    for(int i = 0; i < 9; i++) {
        cout << test3[i] << '\n';
    }*/

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