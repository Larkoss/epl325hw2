#!/bin/bash
thread_num=1
for thread_num in {1..40}
do
	gcc -fopenmp -Werror -Wall -O3 -lm n-body_std.c
    ./a.out 10000 200 $thread_num >> hw2_o3_stats.txt
	gcc -fopenmp -Werror -Wall -O0 -lm n-body_std.c
    ./a.out 10000 200 $thread_num >> hw2_o0_stats.txt
done
