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
        self.mu = None
        self.R = None
        self.k = None
        self.x = None
        self.cp = None
        self._loss = None
        self._ratio_axial_c = None
        self._power = None
        self._pr_ts = None
        self.es = None

        self.def_inlet_radius_ratio = 0.9
        self._d_in = None
        self._d_ex = None
        self._exit_angle = None
        self.exit_vane = None

        self._velocity_diagram = None
        self._max_stage = None
        self._min_stage = None

        self.u_1 = None
        self.u_n = None

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
        self.mu = mu * 1.48817 if unit == US else mu
        self.R = R * 5.38 if unit == US else mu
        self.k = k
        self.x = k / (k - 1)
        self.cp = self.x * self.R

    def set_diameters(self, d_in, d_ex, unit):
        self._d_in = d_in * 25.4e-3 if unit == US else d_in
        self._d_ex = d_ex * 25.4e-3 if unit == US else d_ex

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
        self._exit_angle = alpha * math.pi / 180
        self.exit_vane = False

    def set_exit_vane(self):
        self._exit_angle = None
        self.exit_vane = True

    def set_velocity_diagram(self, diagram):
        self._velocity_diagram = diagram

    def set_number_stage(self, n_min, n_max):
        self._max_stage = n_max
        self._min_stage = n_min

    def get_dh_id(self):
        return self.cp * self._tt * (1 - (1 / self._pr_ts) ** (1 / self.x))

    def get_sum_sq_u(self, n_stage):
        self.u_1 = self._rot * self._d_in
        self.u_n = self._rot * self._d_ex
        if n_stage == 1:
            return self.u_n ** 2
        else:
            sum_sq_u = 0.0
            for i in range(1, n_stage + 1):
                u = (self.u_n - self.u_1) / (n_stage - 1) * (i - 1) + self.u_1
                sum_sq_u += u ** 2
            return sum_sq_u

    def get_q_q2_fls(self, lam):
        q = 0
        q2 = 0
        fls = 0
        if self._velocity_diagram == SYMMETRY_VELOCITY_DIAGRAM:
            q = (lam + 1) / 2
            q2 = (lam - 1) / 2
            fls = 2 - lam

        if self._velocity_diagram == ZERO_EXIT_VELOCITY_DIAGRAM:
            q = 1
            q2 = 0
            fls = 1

        if self._velocity_diagram == IMPULSE_DIAGRAMS:
            q = lam + 0.5
            q2 = lam - 0.5
            fls = 2 * (1 - lam) if lam < 0.5 else 1

        if self._velocity_diagram == ZERO_EXIT_IMPULSE_VELOCITY_DIAGRAM:
            q = 1 if lam >= 0.5 else lam + 0.5
            q2 = 0 if lam >= 0.5 else lam - 0.5
            fls = 1 if lam >= 0.5 else 2 * (1 - lam)
        return q, q2, fls

    def run(self):
        self.es = 0.8
        dh_id = self.get_dh_id()
        p_ex = self._pt / self._pr_ts

        del_ht = self.es * dh_id
        tt_ex = self._tt - del_ht / self.cp

        for n in range(self._min_stage, self._max_stage):
            sum_sq_u = self.get_sum_sq_u(n)
            lam = sum_sq_u / del_ht
            q, q2, fls = self.get_q_q2_fls(lam)
            dv_un = self.u_n / lam
            vu_1n = q * dv_un
            vu_2n = q2 * dv_un
            cot = 1/math.tan(self._exit_angle)
            c1 = (1 + 2 * cot ** 2) * q ** 2
            ci = c1 + (q - 1) ** 2
            d = 2 * cot ** 2 * q ** 2 + (q + lam) ** 2 + (q - lam - 1) ** 2

            re = None
            if n == 1:
                re = self._mf / self.mu / self._d_ex * 2
            else:
                re = self._mf / self.mu / self._d_in * 2

            a = self._loss / cot / re ** 0.2
            a1 = a * (c1 + 2 * d)
            ai = a * (fls * ci + 2 * d)
            b = self._ratio_axial_c * cot ** 2 * q ** 2 + q2 ** 2
            et1 = lam / (lam + a1/2)
            es1 = lam / (lam + (a1 + b)/2)
            eti = lam / (lam + ai / 2)
            esi = lam / (lam + (ai + b) / 2)

            if self.exit_vane:
                cl = 2 * cot ** 2 * q ** 2 + q2 ** 2
                fls_l = 1
                a1l = a1 + a * cl * fls_l
                ail = ai + a * cl * fls_l
                bl = self._ratio_axial_c * cot ** 2 * q ** 2
                et1_l = lam / (lam + a1l / 2)
                es1_l = lam / (lam + (a1l + bl) / 2)
                eti_l = lam / (lam + ail / 2)
                esi_l = lam / (lam + (ail + bl) / 2)
                if n == 1:
                    self.es = es1_l
                else:
                    self.es = 1 / (self.u_1 ** 2 / et1 / sum_sq_u + (1 - (self.u_1 ** 2 + self.u_n ** 2) / sum_sq_u) /
                                   eti_l + self.u_n ** 2 / esi_l / sum_sq_u)
            else:
                if n == 1:
                    self.es = es1
                else:
                    self.es = 1 / (self.u_1 ** 2 / et1 / sum_sq_u +
                                   (1 - (self.u_1 ** 2 + self.u_n ** 2) / sum_sq_u) / eti_l) + \
                              self.u_n ** 2 / esi / sum_sq_u







turbine = PrelimDesign()
turbine.set_stagnation_inlet_parameters(113.1, 2660, US)
turbine.set_regime(53.5, 11400, 0.35, 1.2, US)
turbine.set_fluid(0.376e-4, 53.37, 1.302, US)
turbine.set_power(12900, US)
turbine.set_mean_diameters(22, 24, US)
turbine.set_exit_angle(65)
turbine.set_number_stage(1, 2)
turbine.set_velocity_diagram(SYMMETRY_VELOCITY_DIAGRAM)



