#include<stdio.h>
#include<stdlib.h>
#include<math.h> 


// For brevity or less cumbersome notation we use a 
// macro, now passing functions seems more intuitive
#define FUNC(param1,param2) (*func)(param1,param2)

typedef struct vector{
    int dim;
    double *elements;
} Vector;


Vector *createVector(int dim, double *elements) {
    Vector *vector = malloc(sizeof(Vector));
    vector->dim = dim;
    vector->elements = elements;
    return vector;
}

Vector **createVectorArray(int n, int dim) {
    Vector **arr = malloc(n * sizeof(Vector *));
    int i;
    for (i = 0; i < n; i++) {
        arr[i] = createVector(dim, uniformSeq(dim));
    }
    return arr;
}


double nthHalton(int index, int base) {
            
        double a = 1;
        double val = 0;
        while(index > 0) {
            a = a/base;
            val = val + a*(index % base);
            index = (int) index / base;
        }
        return val;
}
// Make sure to free!
double *haltonSeq(int N, int base) {
    double *sequence = malloc(sizeof(double)*N);
    int i;
    for(i = 0; i < N; i++) {
        sequence[i] = nthHalton(i + 1,base);
    }
    return sequence;
}

// Make sure to free!
double *uniformSeq(int N) {
    double *sequence = malloc(sizeof(double)*N);
    double val;
    int i;
    for(i = 0; i < N; i++) {
        val = (double)rand()/(double)RAND_MAX;
        sequence[i] = val;
    }
    return sequence;
}

void print(double *array, int N) {
    int i;
    for(i = 0; i < N; i++) {
        printf("%f ",array[i]);
    }
}

int main(void) {
    double *seq = haltonSeq(10,2);
    print(seq,10);

    return 0;
}

double myfunction(double *arra) {
    return cos(arra[0])*sin(arra[1]);
}

void linearTransform(Vector **vectorArr, int N, double *lower, double *upper) {
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < vectorArr[i][0]->dim; j++) {
            vectorArr[i][j]->elements[j] = lower[j] + (upper[j] - lower[j]) * vectorArr[i][j]->elements[j];
        }
    }
}
// free ret after use
double * uniformMC(double FUNC(int *points,int dim), int dim, int N, double *lower, double *upper) {
        Vector **points = createVectorArray(N, dim)
        linearTransform(points,lower, upper);
        
        int i;
        double sum = 0;
        for(i = 0; i < N; i++) {
            sum = sum + FUNC(points[i]->elements,dim);
        }
        double factor = 1;
        for(i = 0; i < dim; i++) {
            factor = factor*(upper[i] - lower[i]);
        }
        double variance = 0;
        double mean = factor * sum / N;
        for(i = 0; i < dim; i++) {
            variance = variance  + pow((FUNC(points[i]->elements,dim) - mean),2);
        }
        variance = variance / (N - 1);
        sd = sqrt(variance / N);
        
        freeVectorArr(points);
        
        double * ret = malloc(2 * sizeof(double));
        ret[0] = mean;
        ret[1] = sd;
        return ret;
}


double * miser(void * func, int dim, int N, double * lower, double * upper) {
    int minPoints = N / 10000; // TODO: find good fraction of points
    return miserRecur(func, dim, N, minPoints, lower, upper);
}

double * miserRecur(void * func, int dim, int N, int minPoints, double * lower, double * upper) {
    if (N <= minPoints) {
        // do uniform Monte Carlo Integration on N points
        return uniformMC(func, dim, N, lower, upper);
    } else {
        // recur
        int sampleN = N / 10; // sample 10% of N points in subvolume
        // sample within current bounds of subvolume
        Vector **samplePoints = createVectorArray(sampleN, dim);
        linearTransform(samplePoints, sampleN, lower, upper);
        
        double minSD, sd1;
        int d;
        double * upper1;
        double * lower2;
        int i;
        for (i = 0; i < dim; i++) {
            // bisect along dim[i]
            // compute standard deviations in both bisections
            upper1 = copyArr(upper); // TODO: not an actual function to copy arrays
            lower2 = copyArr(lower); // TODO: not an actual function to copy arrays
            upper1[i] = upper[i] / 2;
            lower2[i] = upper1[i];
            double * bisect1 = uniformMC(func, dim, sampleN, lower, upper1);
            double * bisect2 = uniformMC(func, dim, sampleN, lower2, upper);
            
            if (i == 0) {
                minSD = bisect[1] + bisect[2];
                sd1 = bisect1[1];
                d = 0;
            } else if (bisect1[1] + bisect2[1] < minSD) {
                minSD = bisect1[1] + bisect2[1];
                sd1 = bisect1[1];
                d = i;
            }
            
            // TODO: free upper1, lower2, bisect1, bisect2
        }
        
        // bisect along d
        upper1; // TODO: recopy upper or figure out way to pass d argument recursively
        lower2; // TODO: recopy lower
        upper1[d] = upper[d] / 2;
        lower2[d] = upper1[d];
        int Na = sd1 / minSD * N;
        int Nb = N - Na;
        double * recur1 = miserRecur(func, dim, Na, minPoints, lower, upper1);
        double * recur2 = miserRecur(func, dim, Nb, minPoints, lower2, upper);
        
        // TODO: free upper1 and lower2
        
        double * ret = recur1;
        ret[0] = ret[0] + recur2[0];
        ret[1] = ret[1] + recur2[1];
        return ret;
    }
}