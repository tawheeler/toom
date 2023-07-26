#ifndef PROFILER_H_INCLUDED
#define PROFILER_H_INCLUDED

#include <stdio.h>

#include "platform_metrics.h"
#include "typedefs.h"

#define ArrayCount(Array) (sizeof(Array) / sizeof((Array)[0]))

#ifndef PROFILER
#define PROFILER 1
#endif

#if PROFILER

struct ProfileAnchor
{
    u64 tsc_elapsed_exclusive; // Does _not_ include children
    u64 tsc_elapsed_inclusive; // Does include children
    u64 hit_count;
    char const *label;
};

static struct ProfileAnchor kGlobalProfilerAnchors[4096];
static u32 kGlobalProfilerParent;

struct ProfileBlock
{
    char const *label;
    u64 tsc_start;
    u64 old_tsc_elapsed_inclusive;
    u32 parent_index;
    u32 anchor_index;
};

void InitProfileBlock(struct ProfileBlock *profile_block, char const *label_, u32 anchor_index_)
{
    profile_block->parent_index = kGlobalProfilerParent;
    profile_block->anchor_index = anchor_index_;
    profile_block->label = label_;

    struct ProfileAnchor *anchor = kGlobalProfilerAnchors + anchor_index_;
    profile_block->old_tsc_elapsed_inclusive = anchor->tsc_elapsed_inclusive;

    kGlobalProfilerParent = profile_block->anchor_index;
    profile_block->tsc_start = ReadCPUTimer();
}

void EndProfileBlock(struct ProfileBlock *profile_block)
{
    u64 elapsed = ReadCPUTimer() - profile_block->tsc_start;
    kGlobalProfilerParent = profile_block->parent_index;

    struct ProfileAnchor *parent = kGlobalProfilerAnchors + profile_block->parent_index;
    struct ProfileAnchor *anchor = kGlobalProfilerAnchors + profile_block->anchor_index;

    parent->tsc_elapsed_exclusive -= elapsed;
    anchor->tsc_elapsed_exclusive += elapsed;
    anchor->tsc_elapsed_inclusive = profile_block->old_tsc_elapsed_inclusive + elapsed;
    anchor->hit_count += 1;

    /* NOTE(casey): This write happens every time solely because there is no
        straightforward way in C++ to have the same ease-of-use. In a better programming
        language, it would be simple to have the anchor points gathered and labeled at compile
        time, and this repetative write would be eliminated. */
    anchor->label = profile_block->label;
}

#define NameConcat2(A, B) A##B
#define NameConcat(A, B) NameConcat2(A, B)
#define StartTimeBlock(Name)  \
    struct ProfileBlock Name; \
    InitProfileBlock(Name, __COUNTER__ + 1);
#define EndTimeBlock(Name)
#define ProfilerEndOfCompilationUnit                                \
    static_assert(__COUNTER__ < ArrayCount(kGlobalProfilerAnchors), \
                  "Number of profile points exceeds size of ProfileAnchors array")

static void PrintTimeElapsed(u64 total_tsc_elapsed, struct ProfileAnchor *anchor)
{
    f64 percent = 100.0 * ((f64)anchor->tsc_elapsed_exclusive / (f64)total_tsc_elapsed);
    printf("  %-20s[%6lu]: %10lu (%.2f%%", anchor->label, anchor->hit_count,
           anchor->tsc_elapsed_exclusive, percent);
    if (anchor->tsc_elapsed_inclusive != anchor->tsc_elapsed_exclusive)
    {
        f64 percent_with_children =
            100.0 * ((f64)anchor->tsc_elapsed_inclusive / (f64)total_tsc_elapsed);
        printf(", %.2f%% w/children", percent_with_children);
    }
    printf(")\n");
}

static void PrintAnchorData(u64 total_cpu_elapsed)
{
    for (u32 i_anchor = 0; i_anchor < ArrayCount(kGlobalProfilerAnchors); ++i_anchor)
    {
        struct ProfileAnchor *anchor = kGlobalProfilerAnchors + i_anchor;
        if (anchor->tsc_elapsed_inclusive)
        {
            PrintTimeElapsed(total_cpu_elapsed, anchor);
        }
    }
}

#else

#define TimeBlock(...)
#define PrintAnchorData(...)
#define ProfilerEndOfCompilationUnit

#endif

struct Profiler
{
    u64 tsc_start;
    u64 tsc_end;
};
static struct Profiler kGlobalProfiler;

#define TimeFunction TimeBlock(__func__)

static void BeginProfile(void)
{
    kGlobalProfiler.tsc_start = ReadCPUTimer();
}

static void EndAndPrintProfile()
{
    kGlobalProfiler.tsc_end = ReadCPUTimer();
    u64 cpu_freq = EstimateCPUTimerFreq(100);

    u64 total_cpu_elapsed = kGlobalProfiler.tsc_end - kGlobalProfiler.tsc_start;

    if (cpu_freq)
    {
        printf("\nTotal time: %0.4fms (CPU freq %lu)\n",
               1000.0 * (f64)total_cpu_elapsed / (f64)cpu_freq, cpu_freq);
    }

    PrintAnchorData(total_cpu_elapsed);
}

#endif