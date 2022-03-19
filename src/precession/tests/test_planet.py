from precession.planet import Planet
from precession.config import Config
import yaml


def test_planet_load_dict(snapshot):
    planet = Planet.load(Config(), {
        "M": 30,
    })

    snapshot.assert_match(yaml.safe_dump(
        planet.get_dict()), 'planet_object.yml')


def test_planet_serialize(tmp_path, snapshot):
    planet = Planet.load(Config(), {
        "M": 30,
    })
    planet.save(tmp_path/"planet.yml")
    planet1 = Planet.load(Config(), tmp_path/"planet.yml")

    snapshot.assert_match(yaml.safe_dump(
        planet1.get_dict()), 'planet_object.yml')


def test_planet_fixups(snapshot):
    planet = Planet.load(Config(), {
        "M": 30.,
        "e": 0.1,
        "a": 0.2,
        "T": 0.5
    })
    snapshot.assert_match(f"{planet}", 'planet_string.yml')
    snapshot.assert_match(yaml.safe_dump(
        planet.get_dict()), 'planet_object.yml')
