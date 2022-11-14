// Smart print debugger code
#ifdef __cplusplus
  #include <math.h>
  #include <algorithm>
  #define PRINTLIMIT 1000
  #define print(x);        printScl(__LINE__,#x,x);mexEvalString("drawnow; pause(.005);");mexEvalString("drawnow; pause(.005);");mexEvalString("drawnow; pause(.005);");
  #define printArray(x,N); printArr(__LINE__,#x,x,(int)N);mexEvalString("drawnow; pause(.005);");mexEvalString("drawnow; pause(.005);");mexEvalString("drawnow; pause(.005);");

  // String output
  void printScl(int line, char const *name, char const *       x) {printf("%d string %s:    %s\n",line,name,x);}

  // Scalar outputs
  void printScl(int line, char const *name, bool               x) {printf("%d bool %s:    %s",line,name,x?"true\n":"false\n");}

  void printScl(int line, char const *name, signed char        x) {printf("%d signed char %s:    %hhi\n",line,name,x);}
  void printScl(int line, char const *name, signed short       x) {printf("%d signed short %s:    %hi\n" ,line,name,x);}
  void printScl(int line, char const *name, signed int         x) {printf("%d signed int %s:    %i\n"  ,line,name,x);}
  void printScl(int line, char const *name, signed long        x) {printf("%d signed long %s:    %li\n" ,line,name,x);}
  void printScl(int line, char const *name, signed long long   x) {printf("%d signed long long %s:    %lli\n",line,name,x);}

  void printScl(int line, char const *name, unsigned char      x) {printf("%d unsigned char %s:    %hhu\n",line,name,x);}
  void printScl(int line, char const *name, unsigned short     x) {printf("%d unsigned short %s:    %hu\n" ,line,name,x);}
  void printScl(int line, char const *name, unsigned int       x) {printf("%d unsigned int %s:    %u\n"  ,line,name,x);}
  void printScl(int line, char const *name, unsigned long      x) {printf("%d unsigned long %s:    %lu\n" ,line,name,x);}
  void printScl(int line, char const *name, unsigned long long x) {printf("%d unsigned long long %s:    %llu\n",line,name,x);}

  void printScl(int line, char const *name, float              x) {printf("%d float %s:    %g\n"  ,line,name,x);}
  void printScl(int line, char const *name, double             x) {printf("%d double %s:    %lg\n" ,line,name,x);}
  void printScl(int line, char const *name, long double        x) {printf("%d long double %s:    %Lg\n" ,line,name,x);}

  // Array outputs
  void printArr(int line, char const *name, bool * x, int N) {
    printf("%d bool array %s at %llu:   ",line,name,x);
    bool maxv = false;
    bool minv = true;
    for(int i=0;i<N;i++) {
      if(i<PRINTLIMIT) printf(x[i]?" true":" false",x[i]);
      maxv |= x[i];
      minv &= x[i];
    }
    printf("\nmax: %s, min: %s, length: %d\n",maxv?"true":"false",minv?"true":"false",N);
  }

  // Signed
  void printArr(int line, char const *name, char * x, int N) {
    printf("%d signed char array %s at %llu:   ",line,name,x);
    if(N>0) {
      char maxv = x[0];
      char minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %hhi",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %hhi, min: %hhi, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, signed char * x, int N) {
    printf("%d signed char array %s at %llu:   ",line,name,x);
    if(N>0) {
      signed char maxv = x[0];
      signed char minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %hhi",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %hhi, min: %hhi, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, signed short * x, int N) {
    printf("%d signed short array %s at %llu:   ",line,name,x);
    if(N>0) {
      signed short maxv = x[0];
      signed short minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %hi",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %hi, min: %hi, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, signed int * x, int N) {
    printf("%d signed int array %s at %llu:   ",line,name,x);
    if(N>0) {
      signed int maxv = x[0];
      signed int minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %i",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %i, min: %i, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, signed long * x, int N) {
    printf("%d signed long array %s at %llu:   ",line,name,x);
    if(N>0) {
      signed long maxv = x[0];
      signed long minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %li",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %li, min: %li, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, signed long long * x, int N) {
    printf("%d signed long long array %s at %llu:   ",line,name,x);
    if(N>0) {
      signed long long maxv = x[0];
      signed long long minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %lli",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %lli, min: %lli, length: %d\n",maxv,minv,N);
    }
  }

  // Unsigned
  void printArr(int line, char const *name, unsigned char * x, int N) {
    printf("%d unsigned char array %s at %llu:   ",line,name,x);
    if(N>0) {
      unsigned char maxv = x[0];
      unsigned char minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %hhu",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %hhu, min: %hhu, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, unsigned short * x, int N) {
    printf("%d unsigned short array %s at %llu:   ",line,name,x);
    if(N>0) {
      unsigned short maxv = x[0];
      unsigned short minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %hu",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %hu, min: %hu, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, unsigned int * x, int N) {
    printf("%d unsigned int array %s at %llu:   ",line,name,x);
    if(N>0) {
      unsigned int maxv = x[0];
      unsigned int minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %u",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %u, min: %u, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, unsigned long * x, int N) {
    printf("%d unsigned long array %s at %llu:   ",line,name,x);
    if(N>0) {
      unsigned long maxv = x[0];
      unsigned long minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %lu",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %lu, min: %lu, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, unsigned long long * x, int N) {
    printf("%d unsigned long long array %s at %llu:   ",line,name,x);
    if(N>0) {
      unsigned long long maxv = x[0];
      unsigned long long minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %llu",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %llu, min: %llu, length: %d\n",maxv,minv,N);
    }
  }

  // Floating point precision
  void printArr(int line, char const *name, float * x, int N) {
    printf("%d float array %s at %llu:   ",line,name,x);
    if(N>0) {
      float maxv = x[0];
      float minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %g",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %g, min: %g, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, double * x, int N) {
    printf("%d double array %s at %llu:   ",line,name,x);
    if(N>0) {
      double maxv = x[0];
      double minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %lg",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %lg, min: %lg, length: %d\n",maxv,minv,N);
    }
  }
  void printArr(int line, char const *name, long double * x, int N) {
    printf("%d long double array %s at %llu:   ",line,name,x);
    if(N>0) {
      long double maxv = x[0];
      long double minv = x[0];
      for(int i=0;i<N;i++) {
        if(i<PRINTLIMIT) printf(" %Lg",x[i]);
        maxv = std::max(maxv,x[i]);
        minv = std::min(minv,x[i]);
      }
      printf("\nmax: %Lg, min: %Lg, length: %d\n",maxv,minv,N);
    }
  }
#else
  #define print(x)
  #define printArray(x,N)
#endif

