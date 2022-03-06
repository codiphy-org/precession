from precession.planet import Planet


def test_planet_load_dict():
    planet = Planet.load({
        "mass": 30,
        "radius": 50
    })

    assert planet.mass == 30
    assert planet.radius == 50


def test_planet_serialize(tmp_path):
    planet = Planet.load({
        "mass": 30,
        "radius": 50
    })
    planet.save(tmp_path/"planet.yml")
    planet1 = Planet.load(tmp_path/"planet.yml")

    assert planet1.mass == 30
    assert planet1.radius == 50
