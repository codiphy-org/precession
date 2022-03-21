import math
import numpy as np
import matplotlib.pyplot as plt
import argparse
import logging
import sys
import pandas as pd

from scipy import integrate


logging.basicConfig(level=logging.DEBUG,
                    format='%(levelname)s:%(message)s')


def parse_args():
    parser = argparse.ArgumentParser("Compute precession of orbits")
    parser.add_argument('--config', type=str, default="config.yml",
                        help="yaml config file")
    parser.add_argument('--planet', type=str, default="jupiter.yml",
                        help="planet yml file")
    return parser.parse_args()


def initial_values(config, planet, mercury):
    # initial values are arranged as:
    # [r_p(0), theta_p(0), v_p(0), omega_p(0), r_m(0), theta_m(0), v_m(0), omega_m(0)]
    # note that v = dr_dt and omega = dtheta_dt
    # place the planet and mercury with phase difference of pi
    if (config.onlyMercury):
        return [0., 0., 0., 0., mercury.RMin, 0., 0., (mercury.vMax / mercury.RMin)]
    else:
        return [planet.a, math.pi, 0., (planet.L / planet.a**2), mercury.RMin, 0., 0., (mercury.vMax / mercury.RMin)]


def compute_rhs(t, y, config, planet, mercury):

    # incoming values are arranged as:
    # [r_p(t-), theta_p(t-), v_p(t-), omega_p(t-), r_m(t-), theta_m(t-), v_m(t-), omega_m(t-)]
    # note that v = dr_dt and omega = dtheta_dt
    [r_p, theta_p, v_p, omega_p, r_m, theta_m, v_m, omega_m] = y

    theta_mp = theta_m - theta_p
    cos_theta_mp = math.cos(theta_mp)
    sin_theta_mp = math.sin(theta_mp)
    r_mp = math.sqrt(r_m**2 + r_p**2 - 2 * r_m * r_p * cos_theta_mp)
    r_mp3 = r_mp**3
    alpha_p = planet.GM / r_mp3
    alpha_m = mercury.GM / r_mp3

    d_rp_dt = v_p
    d_thetap_dt = omega_p
    d_vp_dt = (r_p * omega_p**2) - (planet.GMS / r_p**2) - \
        (alpha_m * (r_p - r_m * cos_theta_mp))
    d_omegap_dt = ((-2. * v_p * omega_p) / r_p) - \
        (alpha_m * (r_m / r_p) * sin_theta_mp)

    d_rm_dt = v_m
    d_thetam_dt = omega_m
    d_vm_dt = (r_m * omega_m**2) - (mercury.GMS / r_m**2) - \
        (alpha_p * (r_m - r_p * cos_theta_mp))
    d_omegam_dt = ((-2. * v_m * omega_m) / r_m) + \
        (alpha_p * (r_p / r_m) * sin_theta_mp)

    # return values are arranged as:
    # [d_rp_dt(t), d_thetap_dt(t), d_vp_dt(t), d_omegap_dt(t), d_rm_dt(t), d_thetam_dt(t), d_vm_dt(t), d_omegam_dt(t)]
    return [d_rp_dt, d_thetap_dt, d_vp_dt, d_omegap_dt, d_rm_dt, d_thetam_dt, d_vm_dt, d_omegam_dt]


def solve(config, planet, mercury):
    y = initial_values(config, planet, mercury)
    solution = integrate.solve_ivp(fun=compute_rhs, max_step=config.max_step, t_span=(
        0., config.tEnd), y0=y, atol=config.targetTolerance, args=(config, planet, mercury))
    df = pd.DataFrame(solution.y).transpose()
    df.rename({0: "r_p", 1: "theta_p", 2: "v_p", 3: "omega_p", 4: "r_m",
              5: "theta_m", 6: "v_m", 7: "omega_m"}, axis="columns", inplace=True)
    df["time"] = solution.t
    return df


if __name__ == "__main__":

    from config import Config
    from planet import Planet

    args = parse_args()
    config = Config.load(args.config)

    logging.info(f"{config}")
    planet = Planet.load(config, args.planet)
    logging.info(f"{planet}")
    mercury = Planet.load(config, "mercury.yml")
    logging.info(f"{mercury}")

    df = solve(config, planet, mercury)
