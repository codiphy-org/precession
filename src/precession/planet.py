import yaml
import pathlib
import json
import math


class Planet(object):
    def __init__(self, config):
        self.GMS = 1.
        # mass
        self.M = 0.
        self.name = "unknown"
        # period
        self.T = 1.
        # eccentricity
        self.e = 0.
        # semi major axis
        self.a = 1.
        # configuration
        self.config = config

    def fixup(self):
        self.RMin = self.a * (1 - self.e)
        self.RMax = self.a * (1 + self.e)
        self.vMax = math.sqrt(
            (((1 + self.e) * (1 + self.M)) / self.RMin) * self.GMS)
        self.L = self.a * (1 - self.e) * self.vMax
        if (self.config.decoupledMercury):
            self.GMM = 0.
        else:
            self.GMM = self.GMS * self.M
        if (self.config.grMercury):
            # Constant for General Relativistic correction on Mercury due to Sun
            # This is the value that goes into the force equation term as:
            #  +3*GMS*alpha/r^4
            # And update to potential term in Lagrangian will be something like:
            #  -GMS*alpha/r^3
            # Where 3*alpha value taken from Giordano & Nakanishi is: 1.1e-8
            # self.GMR = self.GMS * (1.1e-8)
            self.GMR = self.GMS / (6.25e+4)**2
        else:
            self.GMR = 0.

    @staticmethod
    def load(config, data):
        if isinstance(data, pathlib.PosixPath):
            data = str(data)
        if isinstance(data, str):
            with open(data, "r") as data_file:
                data = yaml.safe_load(data_file)
        if not isinstance(data, dict):
            raise TypeError(f"data type {type(data)} cannot be loaded")
        planet = Planet(config)
        for k in data:
            setattr(planet, k, data[k])
        planet.fixup()
        return planet

    def get_dict(self):
        data = self.__dict__.copy()
        data.pop("config")
        return data

    def save(self, filename):
        with open(filename, 'w') as file:
            yaml.dump(self.get_dict(), file)

    def __str__(self) -> str:
        return f"planet {self.name} => {', '.join(yaml.safe_dump(self.get_dict()).splitlines())}"
