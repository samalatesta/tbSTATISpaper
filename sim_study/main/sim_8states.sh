#!/bin/bash


# Project name
#$ -P cbs


#Running the main simulation
#Arguments M, p, nstart, initial_iter, N
qsub simsub_8states.qsub 100 .95 10 2000 8
qsub simsub_8states.qsub 250 .95 10 2000 8
qsub simsub_8states.qsub 500 .95 10 2000 8
qsub simsub_8states.qsub 100 .85 10 2000 8
qsub simsub_8states.qsub 250 .85 10 2000 8
qsub simsub_8states.qsub 500 .85 10 2000 8
qsub simsub_8states.qsub 100 .75 10 2000 8
qsub simsub_8states.qsub 250 .75 10 2000 8
qsub simsub_8states.qsub 500 .75 10 2000 8
