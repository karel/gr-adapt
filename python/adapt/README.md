# Tests

To run a specific test use `ctest -R` followed by the test name in the build directory (e.g. `ctest -R qa_lms_filter_ff`).

## Convergence

Convergence tests involve an unknown linear channel whose impulse response is described by the raised cosine<sup>[1](#myfootnote1)</sup>

h = 1/2 * [1 + cos(2pi/W * (n-2))], n = 1,2,3

0, otherwise

where the parameter W controls the amount of amplitude distortion produced by the channel, with the distortion increasing with W. Equivalently, the parameter W controls the eigenvalue spread (i.e., the ratio of the largest eigenvalue to the smallest eigenvalue) of the correlation matrix of the tap inputs of the filter, with the eigenvalue spread increasing with W. Figure 1 presents the learning curves of the filter trained using the LMS algorithm with the step-size parameter μ = 0.0075 and varying W. The filter has M = 11 taps. Figure 2 illustrates the learning curve of LMS at different step-sizes and Figure 3 compares different algorithms.

<img src="./../img/lms_ff_learning_curves_w.svg" alt="Learning curves for the LMS algorithm" title="Figure 1" width="425"><img src="./../img/lms_ff_learning_curves_mu.svg" alt="Learning curves for the LMS algorithm" title="Figure 2" width="425"><img src="./../img/comparison_ff_learning_curves.svg" alt="Learning curves for the different algorithms" title="Figure 3" width="425">

<a name="myfootnote1">1</a>: S. Haykin (2002). Adaptive Filter Theory, 4th Edition, Prentice-Hall.

## Performance

To measure performance, these tests use GNU Radio's [Performance Counters](https://gnuradio.org/doc/doxygen/page_perf_counters.html). By default, GNU Radio will build without Performance Counters enabled. To enable Performance Counters, pass the following flag to cmake:

    -DENABLE_PERFORMANCE_COUNTERS=ON

Given the Performance Counters are enabled in GNU Radio at compile-time, they must also be enabled at run-time. For this, use the GNU Radio preferences file (~/.gnuradio/config.conf) in the section [PerfCounters].

    [PerfCounters]
    on = True
    
### Compiler flags

In addition the the regular optimization level compiler flags, `-ffast-math` can give huge performance boost at the expense of several safety features. This breaks the filters in some cases, be aware.

    add_definitions(-ffast-math)
    add_definitions(-march=native)

### Results

Figures 4 and 5 show performance of the filters provided by gr-adapt with `-ffast-math` disabled and enabled respectively.

<img src="./../img/comparison_xx_performance_i7-6700.svg" alt="Performance of different algorithms" title="Figure 4" width="425"><img src="./../img/comparison_xx_performance_i7-6700_fast-math.svg" alt="Performance of different algorithms with -ffast-math" title="Figure 5" width="425">
