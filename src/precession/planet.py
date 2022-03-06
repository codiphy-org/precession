import yaml
import pathlib


class Planet(object):
    def __init__(self):
        self.mass = 0
        self.name = "unknown"

    @staticmethod
    def load(data):
        if isinstance(data, pathlib.PosixPath):
            data = str(data)
        if isinstance(data, str):
            with open(data, "r") as data_file:
                data = yaml.safe_load(data_file)
        if not isinstance(data, dict):
            raise TypeError(f"data type {type(data)} cannot be loaded")
        planet = Planet()
        for k in data:
            setattr(planet, k, data[k])
        return planet

    def save(self, filename):
        with open(filename, 'w') as file:
            yaml.dump(self.__dict__, file)
