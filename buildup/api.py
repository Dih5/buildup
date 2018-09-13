# -*- coding: utf-8 -*-

"""Python API for buildup"""

import re
import os
from glob import glob

import numpy as np
from scipy import interpolate
from scipy.special import expi
from physdata import xray

_data_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _log_interp_1d(xx, yy, kind='linear'):
    """
    Perform interpolation in log-log scale.
    Args:
        xx (List[float]): x-coordinates of the points.
        yy (List[float]): y-coordinates of the points.
        kind (str or int, optional): The kind of interpolation in the log-log domain. This is passed to
                                     scipy.interpolate.interp1d.
    Returns:
        A function whose call method uses interpolation in log-log scale to find the value at a given point.
    """
    log_x = np.log(xx)
    log_y = np.log(yy)
    # No big difference in efficiency was found when replacing interp1d by UnivariateSpline
    lin_interp = interpolate.interp1d(log_x, log_y, kind=kind, fill_value="extrapolate")
    return lambda zz: np.exp(lin_interp(np.log(zz)))


def get_dose_weight(material="tissue", density=None):
    """
    Get a dose weighting function for the given material, i.e., E*mu_en(E)/rho.

    The mu_en coefficient is fetched from the NIST X-Ray Mass Attenuation Coefficients using the physdata package.

    Args:
        material (str): Name of the material
        density (float): The density of the material. Note the density is not relevant for build-up calculations.

    Returns:
        callable: A function representing the differential fluence to dose conversion coefficient.

    """
    data = np.array(xray.fetch_coefficients(material))
    mu_en = _log_interp_1d(data[:, 0], data[:, 2] * density)
    return lambda energy: energy * mu_en(energy)


class DetectorSingDif:
    """A detector of an energy-differential magnitude from FLUKA"""

    def __init__(self, tab_lis, unit_energy=1.0, unit_integrated=1.0):
        """
        Create a DetectorSingDif instance from a string with its info in the tab.lis format

        Args:
            tab_lis (str): substring of the tab.lis file with the description of results of this detector.
            unit_energy (float): Conversion factor to apply to the FLUKA energies (GeV). E.g., 1E3 to get MeV.
            unit_integrated (float): Conversion factor to apply the scored magnitud. Note the bin size is already
                                     changed by unit_energy.
        """
        if tab_lis == "":
            self.name = ""
            self.number = ""
            return

        r = re.search('# Detector n:\s*([0-9]+)\s*(\S+)', tab_lis)  # Match and parse 1st line
        self.name = r.group(2)
        self.number = r.group(1)

        num_pattern = '-?[0-9]+\.?[0-9]*E[-+]?[0-9]+'
        lines = re.findall((num_pattern + "\s*") * 4, tab_lis)  # Lines with 4 numbers

        # To scale x_min, x_max, y, yrel
        unit_scale = [unit_energy, unit_energy, unit_integrated / unit_energy, 1]
        self.data = []
        for line in lines:
            l_float = map(float, re.findall(num_pattern, line))
            # Scale with the density the magnitudes that depend on it
            self.data.append(list(map(lambda a, b: a * b, l_float, unit_scale)))

        self.data = np.asarray(self.data)

    @classmethod
    def from_data(cls, data):
        """Create a DetectorSingDif from existing data"""
        detector = DetectorSingDif("")
        detector.data = np.asarray(data)
        return detector

    def __repr__(self):
        return "<DetectorSingDif" + str(self.name) + "(" + str(len(self.data)) + ")>"

    def __nonzero__(self):
        """Return truth value (True or False)."""
        return bool(self.data)

    _bool__ = __nonzero__

    def get_xy(self):
        """Return two list with xx and yy pos of the mean points of the histogram"""
        x = []
        y = []
        for l in self.data:
            x.append((l[0] + l[1]) / 2)
            y.append(l[2])
        return x, y

    def get_xerr(self):
        """Return a list with the x interval half-width"""
        xerr = []
        for l in self.data:
            xerr.append((l[1] - l[0]) / 2)
        return xerr

    def get_yerr(self):
        """Return a list with the estimation for the error"""
        yerr = []
        for l in self.data:
            yerr.append(l[2] * l[3] * 0.01)
        return yerr

    def get_norm(self, weight=None, extreme_singular=False):
        """Approximate the value of a norm in the histogram using the mean point"""
        if weight is None:
            return (self.data[:, 2] * (self.data[:, 1] - self.data[:, 0])).sum()
        # weights = np.asarray(list(map(weight, (self.data[:, 0] + self.data[:, 1]) / 2)))
        weights = weight((self.data[:, 0] + self.data[:, 1]) / 2)
        if extreme_singular:
            weights[-1] = weight(self.data[-1, 1])
        return (self.data[:, 2] * weights * (self.data[:, 1] - self.data[:, 0])).sum()

    def get_norm_stat_error(self, weight=None, extreme_singular=False):
        """Approximate the error in a norm due to statistical uncertainty"""
        if weight is None:
            weights = np.ones(len(self.data))
        else:
            weights = weight((self.data[:, 0] + self.data[:, 1]) / 2)
            if extreme_singular:
                weights[-1] = weight(self.data[-1, 1])
        # No need to distinguish singular
        return ((self.data[:, 2] * self.data[:, 3] * 0.01 * weights * (
                self.data[:, 1] - self.data[:, 0])) ** 2).sum() ** 0.5

    def get_norm_hist_error(self, weight=None, extreme_singular=False):
        """Approximate the error in a norm due to binning"""
        if weight is None:
            return 0  # Total counts are exact
        ll = self.data[:-1] if extreme_singular else self.data
        # TODO: If weight is not monotone the max error might not be captured
        weights_mean = weight((ll[:, 1] + ll[:, 0]) / 2)
        weights_low = weight(ll[:, 0])
        weights_hi = weight(ll[:, 1])

        up_est = (ll[:, 2] * np.amax((weights_mean, weights_low, weights_hi), axis=0) * (ll[:, 1] - ll[:, 0])).sum()
        lo_est = (ll[:, 2] * np.amin((weights_mean, weights_low, weights_hi), axis=0) * (ll[:, 1] - ll[:, 0])).sum()
        reg_est = (ll[:, 2] * weights_mean * (ll[:, 1] - ll[:, 0])).sum()
        # Note reg_est is not get_norm, since the singular component is excluded if present
        return max(abs(up_est - reg_est), abs(lo_est - reg_est))

    def get_norm_error(self, weight=None, extreme_singular=False):
        """Get the global error as the 1-norm of the statistical error and the binning bound"""
        return self.get_norm_stat_error(weight=weight, extreme_singular=extreme_singular) + self.get_norm_hist_error(
            weight=weight, extreme_singular=extreme_singular)


def import_single_bdx(tab_lis, unit_energy=1.0, unit_integrated=1.0):
    """
    Return a list of the energy-differential detectors in a tab.lis FLUKA file
    Args:
        tab_lis (str): Path to the tab.lis file.
        unit_energy (float): See DetectorSingDif.
        unit_integrated (float): See DetectorSingDif.

    Yields:
        DetectorSingDif: The next detector in the file.

    """
    detectors_str = re.findall('(# Detector n:.*?)(?=# Detector n:|# double differential distributions|\Z)', tab_lis,
                               re.DOTALL)
    for d in detectors_str:
        yield DetectorSingDif(d, unit_energy=unit_energy, unit_integrated=unit_integrated)


def float_to_str(number):
    """Return a string representation of a float with no rear '.0', as used in the files."""
    s = str(number)
    # No IndexError may rise
    if s[-2:] == '.0':
        s = s[:-2]
    return s


def get_buildup_data(geometry="isotropic", material="lead", weight=None, skip_error=False):
    """Get a build-up coefficient from the npy files"""
    distances = list(np.arange(0.1, 20.001, 0.1))  # Must be the ones in the simulation
    re_pattern = r"buildup-([0-9]+\.?[0-9]*).npy"
    glob_pattern = "buildup-*.npy"
    file_list = glob(os.path.join(_data_path, geometry, material, glob_pattern))
    energies = sorted([float(re.search(re_pattern, f).group(1)) for f in file_list])
    build_up_list = []
    stat_error_list = []
    hist_error_list = []
    if geometry in ["mono", "isotropic"]:
        atten_function = lambda distance: np.exp(-distance)
    elif geometry == "planariso":
        atten_function = lambda distance: -expi(-distance) / 2
    else:
        raise ValueError("Unknown geometry: %s" % geometry)
    for energy in energies:
        weight_energy = 1 if weight is None else weight(energy)
        data_dir = os.path.join(_data_path, geometry, material,
                                "buildup-%s.npy" % float_to_str(energy))
        d_list = [DetectorSingDif.from_data(data) for data in np.load(data_dir)]
        build_up_list.append(
            [d.get_norm(weight, extreme_singular=True) / atten_function(distance) / weight_energy for d, distance in
             zip(d_list, distances)])
        if not skip_error:
            stat_error_list.append(
                [d.get_norm_stat_error(weight, extreme_singular=True) / atten_function(distance) / weight_energy for
                 d, distance in
                 zip(d_list, distances)])
            hist_error_list.append(
                [d.get_norm_hist_error(weight, extreme_singular=True) / atten_function(distance) / weight_energy for
                 d, distance in
                 zip(d_list, distances)])

    if skip_error:
        return BuildUpData(energies, distances, build_up_list)
    else:
        return BuildUpData(energies, distances, build_up_list, stat_error=stat_error_list, hist_error=hist_error_list)


def tab_data_to_bin(geometry="*", material="*", overwrite=False):
    """
    Process tab.lis files in the data directory to create equivalent npy files.

    Write access to the data directory is required for this to work.

    Args:
        geometry (str): Geometry(ies) to consider. Might use a glob pattern.
        material (str):  Material(s) to consider. Might use a glob pattern.
        overwrite (bool):  Whether to overwrite existing files.

    Returns:

    """
    re_pattern = r"buildup-([0-9]+\.?[0-9]*)_26_tab.lis"
    glob_pattern = "buildup-*_26_tab.lis"
    for dir in glob(os.path.join(_data_path, geometry, material)):
        file_list = glob(os.path.join(dir, glob_pattern))
        energies = sorted([float(re.search(re_pattern, f).group(1)) for f in file_list])
        file = 26  # Fixed output unit
        for energy in energies:
            output_file_path = os.path.join(dir, "buildup-%s.npy" % float_to_str(energy))
            if overwrite or not os.path.isfile(output_file_path):
                data_dir = os.path.join(dir, "buildup-%s_%d_tab.lis" % (float_to_str(energy), file,))
                with open(data_dir, 'r') as f:
                    s = f.read()
                np.save(output_file_path,
                        np.asarray(list(d.data for d in import_single_bdx(s, unit_energy=1E3))))


class BuildUpData:
    """
    An object representing a certain build-up factor, which is a function of the distance and the energy and has been
    calculated for a magnitude, a material and a geometry.

    Attributes:
        energies (numpy.ndarray): The energies in the grid.
        distances (numpy.ndarray): The distances in the grid.
        values (numpy.ndarray): The values in the grid.
        stat_error (numpy.ndarray): The error due to statistical uncertainty in the simulations.
        hist_error (numpy.ndarray): The error due to the numerical effect of the grid.

    """

    def __init__(self, energies, distances, values, stat_error=None, hist_error=None):
        self.energies = np.asarray(energies)
        self.distances = np.asarray(distances)
        self.values = np.asarray(values)
        if stat_error is not None:
            self.stat_error = np.asarray(stat_error)
        else:
            self.stat_error = None
        if hist_error is not None:
            self.hist_error = np.asarray(hist_error)
        else:
            self.hist_error = None

    def get_interpolation(self, kx=1, ky=1, error=None):
        """
        Get a function interpolating this build-up data as B(Energy [MeV], distance [mfp]).

        A function extending the build-up factor by means of interpolation will be returned.

        Args:
            kx (int): Spline order in the energy.
            ky (int): Spline order in the distance.
            error (str): Interpolate the error (in some sense) instead of the build-up. Possible values include:
                - "stat": Get the statistical error (the one due to sampling).
                - "hist": Get the bound for the binning error (the one due to making a histogram).
                - "total": Get the sum of both errors.
                - "euclid": Get the euclidean composition of both errors.

        Returns:
            callable: The interpolating B(E, mu*x) function.

        """
        if error is None:
            values = np.asarray(self.values)
        elif error == "stat":
            values = np.asarray(self.stat_error)
        elif error == "hist":
            values = np.asarray(self.hist_error)
        elif error == "total":
            values = np.asarray(self.stat_error) + np.asarray(self.hist_error)
        elif error == "euclid":
            values = np.sqrt(np.asarray(self.stat_error) ** 2 + np.asarray(self.hist_error) ** 2)
        else:
            raise ValueError("Invalidad kwarg error: %s" % str(error))

        return interpolate.RectBivariateSpline(self.energies, self.distances, values, kx=kx, ky=ky)

    def get_interpolated_data(self, energies, distances, kx=1, ky=1):
        """
        Get another BuildUpData instance by changing the mesh by interpolation.

        Args:
            energies (list of float): Energies in the new mesh (in MeV).
            distances (list of float): Distances in the new mesh (in mfp).
            kx (int): Spline order in the energy.
            ky (int): Spline order in the distance.

        Returns:
            BuildUpData: The new instance.

        """
        f = self.get_interpolation(kx=kx, ky=ky)
        e1 = self.get_interpolation(error="stat", kx=kx, ky=ky)
        e2 = self.get_interpolation(error="hist", kx=kx, ky=ky)
        return BuildUpData(energies, distances, f(energies, distances), e1(energies, distances),
                           e2(energies, distances))

    def get_latex_table(self, error=None, transpose=False, number_format="%.2f"):
        """
        Get a LaTeX table describing the data.

        Note you might want to choose a representative mesh with the get_interpolated_data method.

        Args:
            error (str): Whether to display some errors in the table. Possible values include:
                - "stat": Get the statistical error (the one due to sampling).
                - "hist": Get the bound for the binning error (the one due to making a histogram).
                - "total": Get the sum of both errors.
                - "euclid": Get the euclidean composition of both errors.
                " "both": Get both errors independently, first the statistical, then the binning bound.

            transpose (bool): Whether to transpose the table. Default is energy in rows, distance in cols.
            number_format (str): Format string for the numbers.

        Returns:
            str: LaTeX code defining the table.

        """
        # Obtain a matrix with the cells
        if not error:
            # No need for mathematical expression ($~~$) if a single number
            number_format = "%s" % number_format
            values = [[number_format % value for value in row] for row in self.values]
        elif error == "stat":
            number_format = "$%s\\pm%s$" % (number_format, number_format)
            values = [[number_format % value for value in zip(row, error_row)] for row, error_row in
                      zip(self.values, self.stat_error)]
        elif error == "hist":
            number_format = "$%s\\pm%s$" % (number_format, number_format)
            values = [[number_format % value for value in zip(row, error_row)] for row, error_row in
                      zip(self.values, self.stat_error)]
        elif error == "total":
            number_format = "$%s\\pm%s$" % (number_format, number_format)
            error = self.stat_error + self.hist_error
            values = [[number_format % value for value in zip(row, error_row)] for row, error_row in
                      zip(self.values, error)]
        elif error == "euclid":
            number_format = "$%s\\pm%s$" % (number_format, number_format)
            error = np.sqrt(self.stat_error ** 2 + self.hist_error ** 2)
            values = [[number_format % value for value in zip(row, error_row)] for row, error_row in
                      zip(self.values, error)]
        elif error == "both":
            number_format = "$%s\\pm%s\\pm%s$" % (number_format, number_format, number_format)
            values = [[number_format % value for value in zip(*row_tuple)] for row_tuple in
                      zip(self.values, self.stat_error, self.hist_error)]
        else:
            raise ValueError("Invalidad kwarg error: %s" % str(error))

        if transpose:
            values = list(map(list, zip(*values)))
            cols = self.energies
            rows = self.distances
            cols_label = "Energies [MeV]"
            rows_label = "$\\mu x$"
        else:
            cols = self.distances
            rows = self.energies
            cols_label = "Distance [mfp]"
            rows_label = "E [MeV]"

        s = "\\begin{table}\n \\begin{tabular}{%s}\n" % ("c" * (len(cols) + 1))
        s += "& \multicolumn{%d}{c}{%s}\\\\\n" % (len(cols), cols_label)
        s += rows_label + " & " + " & ".join([str(e) for e in cols]) + "\\\\\n"
        for row_label, row in zip(rows, values):
            s += " & ".join([" " + str(row_label)] + row) + "\\\\\n"
        s += "  \\end{tabular}\n \\caption{\\label{tab:my-table}A build-up factor.}\n\\end{table}"
        return s
