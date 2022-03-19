from precession.config import Config
import yaml


def test_config_load_dict():
    config = Config.load({
        "tEnd": 30,
        "maxSteps": 50
    })

    assert config.tEnd == 30
    assert config.maxSteps == 50


def test_config_serialize(tmp_path):
    config = Config.load({
        "tEnd": 30,
        "maxSteps": 50,
    })
    config.save(tmp_path/"config.yml")
    config1 = Config.load(tmp_path/"config.yml")

    assert config1.tEnd == 30
    assert config1.maxSteps == 50


def test_config_fixups(snapshot):
    config = Config.load({
        "tEnd": 30,
        "maxSteps": 50,
    })
    snapshot.assert_match(f"{config}", 'config_string.yml')
    snapshot.assert_match(yaml.safe_dump(config.__dict__), 'config_object.yml')
