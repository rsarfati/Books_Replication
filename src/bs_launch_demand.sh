#!/bin/bash
for B in {1..50}
do
	julia main.jl $B "demand" &
done
