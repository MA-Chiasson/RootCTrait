#!/bin/bash

printf "\nMake structure directory: data, results, results/figures\n"
for i in data results results/figures
do
    mkdir -p "$i"
done