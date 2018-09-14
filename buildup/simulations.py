# -*- coding: utf-8 -*-

"""Create simulations for the buildup package"""

import re
import os
import shutil

import numpy as np
from physdata import xray

from .api import _log_interp_1d

_model_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "model")


def _energy_to_name(energy):
    return str(energy)


def _get_run(template, name, energy, mfp, cycles=20):
    template = re.sub(r"\%\%name\%\%", name, template)
    template = re.sub(r"\%\%energy\%\%", str(energy), template)
    template = re.sub(r"\%\%mfp\%\%", str(mfp), template)
    template = re.sub(r"\%\%cycles\%\%", str(cycles), template)
    return template


# Some usual compounds to be automatically added to inp
_compounds = {}

_compounds["   WATER"] = """* 276 Water liquid H2O
* Chemical Formula:  H -- O -- H
MATERIAL                             1.0                              WATER
COMPOUND         2.0  HYDROGEN       1.0    OXYGEN                    WATER"""

_compounds["CONCRETE"] = """MATERIAL         19.               0.862                              POTASSIU
* Concrete, ordinary as given in NIST x-ray mass attenuation database.
MATERIAL                             2.3                              CONCRETE
COMPOUND     -0.0221  HYDROGEN -0.002484    CARBON  -0.57493    OXYGENCONCRETE
COMPOUND   -0.015208    SODIUM -0.001266  MAGNESIU -0.019953  ALUMINUMCONCRETE
COMPOUND   -0.304627   SILICON -0.010045  POTASSIU -0.042951   CALCIUMCONCRETE
COMPOUND   -0.006435      IRON                                        CONCRETE"""


def create_simulation(path, energy_list, geometry="mono", material_nist=None, material_fluka=None, density=None,
                      energy_to_name=None):
    if not energy_to_name:
        energy_to_name = _energy_to_name
    model_path = os.path.join(_model_path, geometry)

    if material_nist is None:
        raise ValueError("material_nist must be specified")
    if density is None:
        raise ValueError("density must be specified")
    if material_fluka is None:
        raise ValueError("material_fluka must be specified")

    # .flair file
    atten = np.asarray(xray.fetch_coefficients(material_nist, density=density))
    f_atten = _log_interp_1d(atten[:, 0], atten[:, 1])
    with open(os.path.join(model_path, "run.txt"), 'r') as f:
        s = f.read()
    run_list = "\n\n".join((_get_run(s, energy_to_name(e), e, 1.0 / f_atten(e)) for e in energy_list))
    with open(os.path.join(model_path, "buildup.flair"), 'r') as f:
        s = f.read()
    flair_file = re.sub(r"\%\%RUNS\%\%", run_list, s)
    with open(os.path.join(path, "buildup.flair"), 'w') as f:
        f.write(flair_file)

    # .inp file
    with open(os.path.join(model_path, "buildup.inp"), 'r') as f:
        s = f.read()
    material_fluka = material_fluka.upper().rjust(8)  # e.g.: '    LEAD'
    inp_file = re.sub(r"\%\%MATE\%\%", material_fluka, s)
    if material_fluka.upper() in _compounds:
        inp_file = re.sub("GEOEND", "GEOEND\n"+_compounds[material_fluka.upper()], inp_file)
    with open(os.path.join(path, "buildup.inp"), 'w') as f:
        f.write(inp_file)

    shutil.copy(os.path.join(_model_path, "fluscw.f"), path)
