import math

SYMMETRY_VELOCITY_DIAGRAM = 1
ZERO_EXIT_VELOCITY_DIAGRAM = 2
IMPULSE_DIAGRAMS = 3
ZERO_EXIT_IMPULSE_VELOCITY_DIAGRAM = 4

US = 1
SI = 2


class PrelimDesign:
    def __init__(self):

        self._mf = None
        self._rot = None
        self._pt = None
        self._tt = None
        self._mu = None
        self._R = None
        self._k = None
        self._loss = None
        self._ratio_axial_c = None
        self._power = None
        self._pr_ts = None

        self.def_inlet_radius_ratio = 0.9
        self._d_in = None
        self._d_ex = None
        self._exit_angle = None
        self.exit_vane = None

        self._velocity_diagram = None
        self._max_stage = None
        self._min_stage = None

    def set_stagnation_inlet_parameters(self, pt, tt, unit=SI):

        self._pt = pt / 14.504 * 1E+5 if unit == US else pt
        self._tt = tt * 5 / 9 if unit == US else tt

    def set_regime(self, mf, rpm, k_loss, e, unit=SI):
        self._mf = mf / 2.205 if unit == US else mf
        self._rot = rpm * math.pi / 30 if unit == US else rpm
        self._loss = k_loss
        self._ratio_axial_c = e

    def set_power(self, power, unit=SI):
        self._pr_ts = None
        self._power = power * 1.341 * 1000 if unit == US else power

    def set_pr(self, pr_ts):
        self._pr_ts = pr_ts
        self._power = None

    def set_fluid(self, mu, R, k, unit=SI):
        self._mu = mu * 1.48817 if unit == US else mu
        self._R = R * 5.38 if unit == US else mu
        self._k = k

    def set_diameters(self, d_in, d_ex, unit):
        self._d_in = d_in * 25.4 if unit == US else d_in
        self._d_ex = d_ex * 25.4 if unit == US else d_ex

    def set_mean_diameters(self, d_in, d_ex, units=SI):
        self.set_diameters(d_in, d_ex, units)

    def set_tip_diameters(self, d_in, d_ex, exit_radius_ratio, units=SI):
        d_in = (1 + self.def_inlet_radius_ratio) / 2 * d_in
        d_ex = (1 + exit_radius_ratio) / 2 * d_ex
        self.set_diameters(d_in, d_ex, units)

    def set_hub_diameters(self, d_in, d_ex, exit_radius_ratio, units=SI):
        d_in = (1 + self.def_inlet_radius_ratio) / 2 * d_in / self.def_inlet_radius_ratio
        d_ex = (1 + exit_radius_ratio) / 2 * d_ex / exit_radius_ratio
        self.set_diameters(d_in, d_ex, units)

    def set_exit_angle(self, alpha):
        self._exit_angle = alpha
        self.exit_vane = False

    def set_exit_vane(self):
        self._exit_angle = None
        self.exit_vane = True

    def set_velocity_diagram(self, diagram):
        self._velocity_diagram = diagram

    def set_number_stage(self, n_min, n_max):
        self._max_stage = n_max
        self._min_stage = n_min


turbine = PrelimDesign()
turbine.set_stagnation_inlet_parameters(113.1, 2660, US)
turbine.set_regime(53.5, 11400, 0.35, 1.2, US)
turbine.set_fluid(0.376e-4, 53.37, 1.302, US)
turbine.set_power(12900, US)
turbine.set_mean_diameters(22, 24, US)
turbine.set_exit_angle(65)
turbine.set_number_stage(1, 2)
turbine.set_velocity_diagram(SYMMETRY_VELOCITY_DIAGRAM)

#TODO: дописать

