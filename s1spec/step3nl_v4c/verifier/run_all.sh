#!/bin/sh
# ~8 min total with both runs in parallel.  ADD=bcd closes; ADD=abcd does not (see VERIFIER.md).
ADD=bcd TAG=bcd python3 closure.py 1e-3 &
ADD=bcd TAG=bcd python3 closure.py 2e-3 &
wait
ADD="" python3 sanity.py
ADD=abcd python3 cem_probe.py
