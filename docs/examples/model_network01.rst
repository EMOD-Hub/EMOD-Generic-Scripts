===============
model_network01
===============

This example demonstrates the spread of infectivity on a network. Each simulation has a default network of 625 nodes (25-by-25 grid) at 1 of 4 preset levels of network infectivity. Simulations are
implemented as multi-core (four cores per simulation) to demonstrate logging and interaction between nodes hosted on different compute cores (note that cores may or may not be co-located on a common machine.)

An outbreak is initialized by a constant importation pressure of infected individuals. There is no age structure, vital dynamics, or waning. All simulations run for 1000 time steps or until total infectivity falls to zero.

Infected individuals are introduced in one node only (lower-left on the grid). The network infectivity feature implements a gravity-type expression for transmission between nodes. All nodes have equal populations and follow susceptible-infected-recovered dynamics. The value of the exponent in the gravity expression is fixed at 2.0 and the level of network infectivity is varied by setting the coefficient in the gravity expression to different levels.

.. raw:: html
    :file: figures/ref_netinfect_med.svg
