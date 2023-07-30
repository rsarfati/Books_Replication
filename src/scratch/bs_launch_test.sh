#!/bin/bash
for B in {1..4}
do
	julia bs_batch_test.jl $B "demand" &
done
