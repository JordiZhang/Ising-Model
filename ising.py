import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

class ising_model:
    def __init__(self, size, temperature):
        if not isinstance(size, int):
            raise TypeError("Size must be an integer")
        self.size = size
        self.temperature = temperature
        self.rng = np.random.default_rng()

        self.lattice = np.ones((self.size, self.size))
        for i in range(self.size):
            for j in range(self.size):
                if self.rng.uniform(0, 1) > 0.5:
                    self.lattice[i, j] = -1
        self.energy = 0

    def energy_point(self, i, j):
        # calculates energy contribution from a given point
        E = 0
        if i != 0:
            E += -self.lattice[i, j] * self.lattice[i - 1, j]
        if i != self.size - 1:
            E += -self.lattice[i, j] * self.lattice[i + 1, j]
        if j != 0:
            E += -self.lattice[i, j] * self.lattice[i, j - 1]
        if j != self.size - 1:
            E += -self.lattice[i, j] * self.lattice[i, j + 1]
        return E

    def energy_total(self):
        # calculates total energy
        total = 0
        for i in range(self.size):
            for j in range(self.size):
                total += self.energy_point(i, j)
        self.energy = total / 2 # account for double counting

    def glauber_step(self):
        # single iteration of glauber dynamics, includes metropolis algorithm
        for _ in range(self.size*self.size):
            i = self.rng.integers(0, self.size)
            j = self.rng.integers(0, self.size)
            energy_delta = -2 * self.energy_point(i, j)
            #print(energy_delta)
            if energy_delta < 0:
                self.lattice[i, j] = -1 * self.lattice[i, j]
                self.energy += energy_delta
            else:
                r = self.rng.uniform(0, 1)
                p = min(1, np.exp(-energy_delta / self.temperature))
                if r <= p:
                    self.lattice[i, j] = -1 * self.lattice[i, j]
                    self.energy += energy_delta

    def sim_glauber(self):
        # runs a simulation using glauber dynamics and displays it
        fig = plt.figure(figsize=(10, 10))
        im = plt.imshow(self.lattice, interpolation='none', aspect='auto', vmin=-1, vmax=1)
        plt.title("Glauber Dynamics Simulation")

        def animate_func(i):
            self.glauber_step()
            print(self.energy)
            im.set_array(self.lattice)
            return [im]

        anim = animation.FuncAnimation(fig, animate_func, interval=1, blit=True)
        plt.show()

    def kawasaki_step(self):
        for _ in range(int(self.size*self.size/2)):
            i1 = self.rng.integers(0, self.size)
            j1 = self.rng.integers(0, self.size)
            # makes sure the 2 points are distinct from each other
            while True:
                i2 = self.rng.integers(0, self.size)
                j2 = self.rng.integers(0, self.size)
                if i1 != i2 or j1 != j2:
                    break
            S1 = self.lattice[i1, j1]
            S2 = self.lattice[i2, j2]

            # only if they have different spins does swapping make a difference
            if S1 != S2:
                # if nearest neighbours
                if (np.absolute(i1-i2) == 1 and j1 == j2) or (np.absolute(j1-j2) == 1 and i1 == i2):
                    # energy delta from 1 flip, a second flip, + 4
                    energy_delta = -2*self.energy_point(i1, j1) + -2*self.energy_point(i2, j2) + 4
                else:
                    energy_delta = -2*self.energy_point(i1, j1) + -2*self.energy_point(i2, j2)

                # metropolis
                # print(energy_delta)
                if energy_delta < 0:
                    # swaps them
                    self.lattice[i1, j1] = S2
                    self.lattice[i2, j2] = S1
                    self.energy += energy_delta
                else:
                    r = self.rng.uniform(0, 1)
                    p = min(1, np.exp(-energy_delta / self.temperature))
                    if r <= p:
                        self.lattice[i1, j1] = S2
                        self.lattice[i2, j2] = S1
                        self.energy += energy_delta

    def sim_kawasaki(self):
        # runs a simulation using kawasaki dynamics and displays it

        fig = plt.figure(figsize=(10, 10))
        im = plt.imshow(self.lattice, interpolation='none', aspect='auto', vmin=-1, vmax=1)
        plt.title("Kawasaki Dynamics Simulation")
        def animate_func(i):
            self.kawasaki_step()
            print(self.energy)
            im.set_array(self.lattice)
            return [im]

        anim = animation.FuncAnimation(fig, animate_func, interval=1, blit=True)
        plt.show()

    def measurement(self, dynamics, temperature):
        # dynamics: k, kawasaki. g, glauber
        # Faster equilibrium
        if dynamics == "g":
            self.lattice = np.ones((self.size, self.size))
            self.energy_total()
        if dynamics == "k":
            self.lattice = np.ones((self.size, self.size))
            self.lattice[int(self.size/2):-1, :] = self.lattice[int(self.size/2):-1, :]*-1
            self.energy_total()

        for i in range(len(temperature)):
            self.temperature = temperature[i]

            if dynamics == "g":
                print("Glauber Simulation in progress. T = " + str(self.temperature))
                for _ in range(100):
                    self.glauber_step()

            if dynamics == "k":
                print("Kawasaki Simulation in progress. T = " + str(self.temperature))
                for _ in range(100):
                    self.kawasaki_step()

            S = np.sum(self.lattice)
            M = [S]
            M2 = [S*S]
            E = [self.energy]
            E2 = [self.energy*self.energy]

            # takes 1000 measurements and saves them
            for i in range(1000):
                if i % 100 == 0:
                    print("Progress: " + str(i / 10))
                for j in range(15):
                    if dynamics == "g":
                        self.glauber_step()
                    if dynamics == "k":
                        self.kawasaki_step()
                S = np.sum(self.lattice)
                M.append(S)
                M2.append(S*S)
                E.append(self.energy)
                E2.append(self.energy*self.energy)

            if dynamics == "g":
                print("Glauber Simulation Finished. T = " + str(self.temperature))
            if dynamics == "k":
                print("Kawasaki Simulation Finished. T = " + str(self.temperature))

            M = np.array(M).reshape((-1,1))
            M2 = np.array(M2).reshape((-1, 1))
            E = np.array(E).reshape((-1,1))
            E2 = np.array(E2).reshape((-1, 1))

            filename = dynamics + str(self.temperature) + ".csv"
            data = np.concatenate((M, M2, E, E2), axis = 1)
            np.savetxt(fname = filename, X = data, delimiter=",", header = "M, M^2, E, E^2")
