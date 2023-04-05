#ifndef MERGE_SORT_H
#define MERGE_SORT_H

#include <algorithm>
#include <vector>
#include <omp.h>

#ifndef TIE_THRESHOLD
#define TIE_THRESHOLD (1<<14)
#endif

#ifndef SERIAL_CUTOFF
#define SERIAL_CUTOFF (1<<18)
#endif

/**
 * Parallel Merge Sort
 *
 * Taken from https://stackoverflow.com/a/13811342
 */
template <class RandomIT, class Compare>
void mergeSortRecursive(RandomIT left, RandomIT right, Compare comp) 
{
    if (left < right) 
    {
        if (std::distance(left, right) >= SERIAL_CUTOFF) 
	{
            RandomIT mid = left+std::distance(left, right)/2;
            #pragma omp taskgroup
            {
                #pragma omp task untied if(right-left >= TIE_THRESHOLD)
                mergeSortRecursive(left, mid, comp);
                #pragma omp task untied if(right-left >= TIE_THRESHOLD)
                mergeSortRecursive(mid, right, comp);
                #pragma omp taskyield
            }
            std::inplace_merge(left, mid, right, comp);
        } 
	else 
	{
            std::sort(left, right, comp);
        }
    }
}

template <class RandomIT, class Compare>
void MergeSort(RandomIT begin, RandomIT end, Compare comp) 
{
    #pragma omp parallel
    #pragma omp single
    mergeSortRecursive(begin, end, comp);
}

#endif