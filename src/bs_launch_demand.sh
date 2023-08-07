#!/bin/bash
for B in {1..20}
do
	julia main.jl $B "demand" &
done
