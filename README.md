# Ising Model 
---
This project implements a 2D Ising model, simulating ferromagnetic spin systems on a lattice using both Glauber and Kawasaki dynamics.

## Installing
Download ising.py and import.
```
from ising import IsingModel
```

## Usage
Create an IsingModel object as needed.
```
# ising model parameters
size = 50
temperature = 0.1

model = IsingModel(size, temperature)
```
Initialize the total energy of the system.
```
model.energy_total()
```
Simulate using Glauber or Kawasaki Dynamics. The simulation is automatically animated on matplotlib.
```
model.sim_glauber()
```
or
```
model.sim_kawasaki()
```
