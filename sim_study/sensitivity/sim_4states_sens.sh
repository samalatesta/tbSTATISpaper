#!/bin/bash


# Project name
#$ -P cbs

#join error and output stream
#$ -j y

#Running the main simulation
#Arguments M, p, nstart, initial_iter, N
qsub simsub_4states_sens.qsub 100 .95 10 2000 4
qsub simsub_4states_sens.qsub 250 .95 10 2000 4
qsub simsub_4states_sens.qsub 500 .95 10 2000 4
qsub simsub_4states_sens.qsub 100 .85 10 2000 4
qsub simsub_4states_sens.qsub 250 .85 10 2000 4
qsub simsub_4states_sens.qsub 500 .85 10 2000 4
qsub simsub_4states_sens.qsub 100 .75 10 2000 4
qsub simsub_4states_sens.qsub 250 .75 10 2000 4
qsub simsub_4states_sens.qsub 500 .75 10 2000 4
