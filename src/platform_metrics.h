#ifndef PLATFORM_METRICS_C_INCLUDED
#define PLATFORM_METRICS_C_INCLUDED

// This platform metrics file is based on the performance aware programming series by Casey Muratori.

#include "typedefs.h"

#include <x86intrin.h>

// NOTE: Linux does not provide a way to obtain the OS timer frequency.
u64 GetOSTimerFreq();

u64 ReadOSTimer();

u64 ReadCPUTimer();

f64 GetElapsedCPUTimeMs(u64 timer_start, u64 timer_end, u64 cpu_timer_freq);

u64 EstimateCPUTimerFreq(u64 wait_time_ms);

#endif