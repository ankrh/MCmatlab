# Enabling OpenMP on macOS for MCmatlab

Thanks to Will Grissom who figured this method out for his tool kpTx (https://github.com/wgrissom)

The following strategy has been tested on:
- MATLAB R2022a running on macOS Ventura 13.0
- MATLAB R2023b running on macOS Sonoma 14.0
- MATLAB R2024b running on macOS Sequoia 15.0

## Prerequisites

1. Install XCode from the App Store and subsequently install the Apple Command Line Tools:
   (enter the following command in the terminal)
   ```
   xcode-select --install
   ```

2. Install Homebrew:
   (enter the following command in the terminal)
   ```
   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
   ```

3. Use Homebrew to install LLVM and OpenMP:
   (enter the following command in the terminal)
   ```
   brew install llvm libomp
   ```

## Compilation Steps

1. Configure MATLAB to use the custom XML configuration file:
   (enter the following command on the MATLAB command line)
   ```
   mex -setup:+MCmatlab/src/clang_openmp_maci64.xml C
   ```
   (clang_openmp_maci64.xml is included in this repository)

2. Compile the MEX file with OpenMP support:
   (enter the following command on the MATLAB command line)
   ```
   mex COPTIMFLAGS='$COPTIMFLAGS -O3 -ffast-math -fno-finite-math-only -fopenmp -std=c11 -Wall' LDOPTIMFLAGS='$LDOPTIMFLAGS -O3 -ffast-math -fno-finite-math-only -fopenmp -std=c11 -Wall' -outdir +MCmatlab/@model/private ./+MCmatlab/src/MCmatlab.c
   ```

## Optimization Flags Explanation

The command line flags override the XML configuration and provide:
- `-O3`: Highest standard optimization level
- `-ffast-math`: Aggressive floating-point optimizations
- `-fno-finite-math-only`: Don't assume arguments and results are not NaN/Inf (which -ffast-math would otherwise do)
- `-fopenmp`: Enable OpenMP multithreading
- `-std=c11`: Use C11 language standard
- `-Wall`: Enable all common warnings

The following sites were helpful in figuring this out:

https://github.com/wgrissom/kpTx, and therein
https://stackoverflow.com/questions/37362414/openmp-with-mex-in-matlab-on-mac
https://stackoverflow.com/questions/43555410/enable-openmp-support-in-clang-in-mac-os-x-sierra-mojave