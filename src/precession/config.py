import yaml
import pathlib
import json


class Config(object):
    def __init__(self):
        self.tEnd = 1
        self.maxSteps = 1
        self.batchSteps = 1
        self.filePrefix = "precession"
        self.onlyMercury = False
        self.decoupledMercury = False
        self.grMercury = True
        self.enableConvergenceTest = True
        self.GMS = 2
        self.G = 2

    @staticmethod
    def load(data):
        if isinstance(data, pathlib.PosixPath):
            data = str(data)
        if isinstance(data, str):
            with open(data, "r") as data_file:
                data = yaml.safe_load(data_file)
        if not isinstance(data, dict):
            raise TypeError(f"data type {type(data)} cannot be loaded")
        planet = Config()
        for k in data:
            setattr(planet, k, data[k])
        return planet

    def save(self, filename):
        with open(filename, 'w') as file:
            yaml.dump(self.__dict__, file)

    def __str__(self) -> str:
        return f"config => {', '.join(yaml.safe_dump(self.__dict__).splitlines())}"
