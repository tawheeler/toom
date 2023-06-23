#include "platform_metrics.h"

#include <sys/time.h>

u64 GetOSTimerFreq() { return 1000000; }

u64 ReadOSTimer()
{
    // NOTE(casey): The "struct" keyword is not necessary here when compiling in C++,
    // but just in case anyone is using this file from C, I include it.
    struct timeval value;
    gettimeofday(&value, 0);

    u64 result = GetOSTimerFreq() * (u64)value.tv_sec + (u64)value.tv_usec;
    return result;
}

u64 ReadCPUTimer() { return __rdtsc(); }

f64 GetElapsedCPUTimeMs(u64 timer_start, u64 timer_end, u64 cpu_timer_freq)
{
    u64 timer_delta = timer_end - timer_start;
    f64 elapsed_ms = (timer_delta * 1000.0) / cpu_timer_freq;
    return elapsed_ms;
}

u64 EstimateCPUTimerFreq(u64 wait_time_ms)
{
    u64 os_freq = GetOSTimerFreq();
    u64 os_wait_time = os_freq * wait_time_ms / 1000;

    u64 cpu_start = ReadCPUTimer();
    u64 os_start = ReadOSTimer();
    u64 os_end = 0;
    u64 os_elapsed = 0;
    while (os_elapsed < os_wait_time)
    {
        os_end = ReadOSTimer();
        os_elapsed = os_end - os_start;
    }

    u64 cpu_end = ReadCPUTimer();
    u64 cpu_elapsed = cpu_end - cpu_start;
    u64 cpu_freq = 0;
    if (os_elapsed)
    {
        cpu_freq = os_freq * cpu_elapsed / os_elapsed;
    }
    return cpu_freq;
}