import math

SYMMETRY_VELOCITY_DIAGRAM = 1
ZERO_EXIT_VELOCITY_DIAGRAM = 2
IMPULSE_DIAGRAMS = 3
ZERO_EXIT_IMPULSE_VELOCITY_DIAGRAM = 4

US = 1
SI = 2


class PrelimDesign:
    def __init__(self):

        # Исходные данные
        self.mf = None
        self.rot = None
        self.pt_in = None
        self.tt_in = None
        self.mu = None
        self.R = None
        self.k = None
        self.loss = None
        self.ratio_axial_c = None
        self.power = None
        self.pr_ts = None

        self.d_in = None
        self.d_ex = None
        self.exit_angle = None
        self.exit_vane = None
        self.velocity_diagram = None

        self.max_stage = None
        self.min_stage = None

        # внутренние переменные
        self.x = None
        self.cp = None
        self.p_ex = None
        self.t_ex = None
        self.ro = None
        self.tt_ex = None
        self.es = 0.8

        self.sum_sq_u = None
        self.lam = None
        self.q = None
        self.q2 = None
        self.fls = None
        self.re = None
        self.zz = None

        self.dv_un = None
        self.vu_1n = None
        self.vu_2n = None

        self.vx_nd = None
        self.vx_n = None
        self.v2_n = None
        self.v2_nr = None

        self.inlet_radius_ratio = 0.9
        self.exit_radius_ratio = None
        self._d_in_ = None
        self._d_ex_ = None
        self.func_d = None

        self.u_1 = None
        self.u_n = None

        self.conv = 0.01

    def print_task(self):
        print(5*'\t' + 'TASK INFO' + 5*'\t')
        print(f'FLUID: R = {self.R:.1f}, mu = {self.mu:.3e}, k = {self.k:.2f}')
        print(f'INLET: mf = {self.mf:.1f} кг/с, pt_in = {self.pt_in*1e-5:.2f} бар, tt_in = {self.tt_in:.0f} К')
        task_string = f'REGIME: rpm = {self.rot * 30/math.pi:.0f} об/мин. '
        if self.power is None:
            print(task_string + f'pr_ts = {self.pr_ts}')
        else:
            print(task_string + f'power = {self.power/1e+3:.0f} kВт')
        print(f'DIAMETERS: inlet = {self.d_in * 1e+3:.1f} мм exit = {self.d_ex * 1e+3:.1f} мм')
        if not self.exit_vane:
            print(f'ANGLE: exit_angle = {self.exit_angle * 180/math.pi:.1f}')

    def set_stagnation_inlet_parameters(self, pt, tt, unit=SI):

        self.pt_in = pt * 6894.76 if unit == US else pt
        self.tt_in = tt * 5 / 9 if unit == US else tt

    def set_regime(self, mf, rpm, k_loss, e, unit=SI):
        self.mf = mf / 2.205 if unit == US else mf
        self.rot = rpm * math.pi / 30 if unit == US else rpm
        self.loss = k_loss
        self.ratio_axial_c = e

    def set_power(self, power, unit=SI):
        self.pr_ts = None
        self.power = power * 1.341 * 1000 if unit == US else power

    def set_pr(self, pr_ts):
        self.pr_ts = pr_ts
        self.power = None

    def set_fluid(self, mu, R, k, unit=SI):
        self.mu = mu * 1.48817 if unit == US else mu
        self.R = R * 5.38 if unit == US else mu
        self.k = k
        self.x = k / (k - 1)
        self.cp = self.x * self.R

    def set_diameters(self, d_in, d_ex, unit):
        self.d_in = d_in / 39.37 if unit == US else d_in
        self.d_ex = d_ex / 39.37 if unit == US else d_ex

    def set_mean_diameters(self, d_in, d_ex, units=SI):
        self.set_diameters(d_in, d_ex, units)

    def set_tip_diameters(self, d_in, d_ex, exit_radius_ratio, units=SI):
        if self.func_d is None:
            self.func_d = self.set_tip_diameters

        self._d_in_ = d_in
        self._d_ex_ = d_ex
        self.exit_radius_ratio = exit_radius_ratio
        d_in = (1 + self.inlet_radius_ratio) / 2 * d_in
        d_ex = (1 + exit_radius_ratio) / 2 * d_ex
        self.set_diameters(d_in, d_ex, units)

    def set_hub_diameters(self, d_in, d_ex, exit_radius_ratio, units=SI):
        if self.func_d is None:
            self.func_d = self.set_hub_diameters

        self._d_in_ = d_in
        self._d_ex_ = d_ex
        self.exit_radius_ratio = exit_radius_ratio
        d_in = (1 + self.inlet_radius_ratio) / 2 * d_in / self.inlet_radius_ratio
        d_ex = (1 + exit_radius_ratio) / 2 * d_ex / exit_radius_ratio
        self.set_diameters(d_in, d_ex, units)

    def set_exit_angle(self, alpha):
        self.exit_angle = alpha * math.pi / 180
        self.exit_vane = False

    def set_exit_vane(self):
        self.exit_angle = None
        self.exit_vane = True

    def set_velocity_diagram(self, diagram):
        self.velocity_diagram = diagram

    def set_number_stage(self, n_min, n_max):
        self.max_stage = n_max
        self.min_stage = n_min

    def get_dh_id(self):
        return self.cp * self.tt_in * (1 - (1 / self.pr_ts) ** (1 / self.x))

    def get_del_ht(self):
        return self.power / self.mf

    def get_sum_sq_u(self, n_stage):
        self.u_1 = self.rot * self.d_in
        self.u_n = self.rot * self.d_ex
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
        if self.velocity_diagram == SYMMETRY_VELOCITY_DIAGRAM:
            q = (lam + 1) / 2
            q2 = (lam - 1) / 2
            fls = 2 - lam

        if self.velocity_diagram == ZERO_EXIT_VELOCITY_DIAGRAM:
            q = 1
            q2 = 0
            fls = 1

        if self.velocity_diagram == IMPULSE_DIAGRAMS:
            q = lam + 0.5
            q2 = lam - 0.5
            fls = 2 * (1 - lam) if lam < 0.5 else 1

        if self.velocity_diagram == ZERO_EXIT_IMPULSE_VELOCITY_DIAGRAM:
            q = 1 if lam >= 0.5 else lam + 0.5
            q2 = 0 if lam >= 0.5 else lam - 0.5
            fls = 1 if lam >= 0.5 else 2 * (1 - lam)
        return q, q2, fls

    def run_pr_ts(self, n_stage):
        dh_id = self.get_dh_id()
        del_ht = self.es * dh_id

        _es_ = self.es
        while abs((es := self.get_es(del_ht, n_stage)) - _es_) > 1e-4:
            self.es = (es + _es_) / 2
            _es_ = es
            del_ht = self.es * dh_id

        self.p_ex = self.pt_in / self.pr_ts
        self.tt_ex = self.tt_in - del_ht / self.cp

    def run_power(self, n_stage):
        del_ht = self.get_del_ht()

        self.tt_ex = self.tt_in - del_ht / self.cp
        if self.tt_ex < 0.0:
            raise Exception('Недостаточно энергии на входе!')

        self.es = self.get_es(del_ht, n_stage)
        if zz := del_ht / (self.cp * self.tt_in * self.es) >= 1:
            raise Exception('Недостаточно распологаемого теплоперепада')

        self.p_ex = self.pt_in * (1 - zz) ** self.x

        cot = 1 / math.tan(self.exit_angle)

        self.vx_nd = self.q * cot * self.dv_un
        self.vx_n = self.vx_nd * self.ratio_axial_c ** 0.5
        self.v2_n = (self.vu_2n ** 2 + self.vx_n ** 2) ** 0.5
        self.v2_nr = self.v2_n

        if self.exit_vane:
            self.v2_n = self.vx_n

        self.t_ex = self.tt_ex - self.v2_n ** 2 / (self.mf * self.cp)
        if self.t_ex < 0.0:
            raise Exception('Отрицательная температцра на выходе')

        self.ro = self.p_ex / (self.R * self.t_ex)
        # если задано соотношени радиусов go to 8
        # TODO:
        if self.exit_angle is None:
            area_exit = self.mf / self.ro / self.vx_n
            aa = area_exit / math.pi / self.d_ex ** 2
            if aa >= 1:
                print('Недостаточная площадь на выходе')
                self.inlet_radius_ratio = 0.9
                self.func_d(self._d_in_, self._d_ex_, 0.7)
                return
        rr_exit = (1 - aa) / (1 + aa)
        # TODO: реализация логики для задания диаметров на ступице и переферии

        # if n_stage == 1:
        #     self.conv = 1e-4
        # if abs(rr_exit - self.exit_radius_ratio) < self.conv:
        #     # TODO:
        #     pass








    def get_es(self, del_ht, n_stage):
        self.sum_sq_u = self.get_sum_sq_u(n_stage)
        self.lam = self.sum_sq_u / del_ht
        self.q, self.q2, self.fls = self.get_q_q2_fls(self.lam)
        self.dv_un = self.u_n / self.lam
        self.vu_1n = self.q * self.dv_un
        self.vu_2n = self.q2 * self.dv_un

        cot = 1 / math.tan(self.exit_angle)

        c1 = (1 + 2 * cot ** 2) * self.q ** 2
        ci = c1 + (self.q - 1) ** 2
        d = 2 * cot ** 2 * self.q ** 2 + (self.q + self.lam) ** 2 + (self.q - self.lam - 1) ** 2
        self.re = None
        if n_stage == 1:
            self.re = self.mf / self.mu / self.d_ex * 2
        else:
            self.re = self.mf / self.mu / self.d_in * 2
        a = self.loss / cot / self.re ** 0.2
        a1 = a * (c1 + 2 * d)
        ai = a * (self.fls * ci + 2 * d)
        b = self.ratio_axial_c * cot ** 2 * self.q ** 2 + self.q2 ** 2
        et1 = self.lam / (self.lam + a1 / 2)
        es1 = self.lam / (self.lam + (a1 + b) / 2)
        eti = self.lam / (self.lam + ai / 2)
        esi = self.lam / (self.lam + (ai + b) / 2)
        if self.exit_vane:
            cl = 2 * cot ** 2 * self.q ** 2 + self.q2 ** 2
            fls_l = 1
            a1l = a1 + a * cl * fls_l
            ail = ai + a * cl * fls_l
            bl = self.ratio_axial_c * cot ** 2 * self.q ** 2
            et1_l = self.lam / (self.lam + a1l / 2)
            es1_l = self.lam / (self.lam + (a1l + bl) / 2)
            eti_l = self.lam / (self.lam + ail / 2)
            esi_l = self.lam / (self.lam + (ail + bl) / 2)
            if n_stage == 1:
                return es1_l
            else:
                return 1 / (self.u_1 ** 2 / et1 / self.sum_sq_u + (1 - (self.u_1 ** 2 + self.u_n ** 2) /
                            self.sum_sq_u) / eti_l + self.u_n ** 2 / esi_l / self.sum_sq_u)
        else:
            if n_stage == 1:
                return es1
            else:
                return 1 / (self.u_1 ** 2 / et1 / self.sum_sq_u + (1 - (self.u_1 ** 2 + self.u_n ** 2) /
                            self.sum_sq_u) / eti) + self.u_n ** 2 / esi / self.sum_sq_u


turbine = PrelimDesign()
turbine.set_stagnation_inlet_parameters(pt=113.1, tt=2660, unit=US)
turbine.set_regime(mf=53.5, rpm=11400, k_loss=0.35, e=1.2, unit=US)
turbine.set_fluid(mu=0.376e-4, R=53.37, k=1.302, unit=US)
# turbine.set_power(power=12900, unit=US)
turbine.set_pr(pr_ts=3.66)
turbine.set_mean_diameters(d_in=22, d_ex=24, units=US)
turbine.set_exit_angle(alpha=65)
turbine.set_number_stage(n_min=1, n_max=2)
turbine.set_velocity_diagram(SYMMETRY_VELOCITY_DIAGRAM)
# turbine.run_power(1)
turbine.run_pr_ts(1)
