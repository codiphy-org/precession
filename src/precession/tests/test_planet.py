from precession.planet import Planet
from precession.config import Config
import yaml


def test_planet_load_dict():
    planet = Planet.load(Config(), {
        "M": 30,
        "R": 50
    })

    assert planet.M == 30
    assert planet.R == 50


def test_planet_serialize(tmp_path):
    planet = Planet.load(Config(), {
        "M": 30,
        "R": 50,
    })
    planet.save(tmp_path/"planet.yml")
    planet1 = Planet.load(Config(), tmp_path/"planet.yml")

    assert planet1.M == 30
    assert planet1.R == 50


def test_planet_fixups(snapshot):
    planet = Planet.load(Config(), {
        "M": 30.,
        "R": 50.,
        "e": 0.1,
        "a": 0.2,
        "T": 0.5
    })
    snapshot.assert_match(f"{planet}", 'planet_string.yml')
    snapshot.assert_match(yaml.safe_dump(
        planet.get_dict()), 'planet_object.yml')
