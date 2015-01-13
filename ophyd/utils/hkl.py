# vi: ts=4 sw=4 sts=4 expandtab
'''
:mod:`ophyd.utils.hkl` - HKL calculation utilities
==================================================

.. module:: ophyd.utils.hkl
   :synopsis:

'''

from __future__ import print_function
import numpy as np

import sys
try:
    from gi.repository import Hkl as HklModule
    from gi.repository import GLib
except ImportError as ex:
    print('[!!] Failed to import Hkl library; diffractometer support'
          ' disabled (%s)' % ex,
          file=sys.stderr)

    HklModule = None

from ..controls import PseudoPositioner


def new_detector(dtype=0):
    '''
    Create a new HKL-library detector
    '''
    return HklModule.Detector.factory_new(HklModule.DetectorType(dtype))


if HklModule:
    DIFF_TYPES = tuple(sorted(HklModule.factories().keys()))
    UserUnits = HklModule.UnitEnum.USER
else:
    DIFF_TYPES = ()


class UsingEngine(object):
    """
    Context manager that uses a calculation engine temporarily (i.e., for the
    duration of the context manager)
    """
    def __init__(self, calc, engine):
        self.calc = calc

    def __enter__(self):
        self.old_engine = self.calc.engine

    def __exit__(self, type_, value, traceback):
        if self.old_engine is not None:
            self.calc.engine = self.old_engine


class HklCalc(object):
    def __init__(self, dtype, engine='hkl',
                 sample='main', lattice=None,
                 degrees=True):
        self._engine = None  # set below with property
        self._detector = new_detector()
        self._degrees = bool(degrees)
        self._sample = None
        self._samples = {}

        try:
            self._factory = HklModule.factories()[dtype]
        except KeyError:
            raise ValueError('Invalid diffractometer type (%s)'
                             'Choose from: %s' % (dtype, ', '.join(DIFF_TYPES)))

        self._geometry = self._factory.create_new_geometry()
        self._engine_list = self._factory.create_new_engine_list()
        self._solutions = None

        if sample is not None:
            self.add_sample(sample, lattice=lattice)

        self.engine = engine

    @property
    def engine(self):
        return self._engine

    @engine.setter
    def engine(self, engine):
        if engine is self._engine:
            return

        if isinstance(engine, HklModule.Engine):
            self._engine = engine
        else:
            engines = self.engines
            try:
                self._engine = engines[engine]
            except KeyError:
                raise ValueError('Unknown engine name or type')

        self._re_init()

    @property
    def sample_name(self):
        '''
        The name of the currently selected sample
        '''
        return self._sample.name_get()

    @sample_name.setter
    def sample_name(self, new_name):
        sample = self._sample
        current = sample.name_get()
        if new_name in self._samples:
            raise ValueError('Sample with that name already exists')

        sample.name_set(new_name)

        del self._samples[current]
        self._samples[new_name] = sample

    @property
    def sample(self):
        return self._sample

    @sample.setter
    def sample(self, sample):
        if sample is self._sample:
            return

        if isinstance(sample, HklModule.Sample):
            if sample not in self._samples.values():
                self.add_sample(sample)
            else:
                self._sample = sample
                self._re_init()
        elif sample in self._samples:
            name = sample
            self._sample = self._samples[name]
            self._re_init()
        else:
            raise ValueError('Unknown sample type (expected Hkl.Sample)')

    def add_sample(self, name, lattice=None, select=True):
        if isinstance(name, HklModule.Sample):
            sample = name
            name = sample.name_get()
        else:
            sample = HklModule.Sample.new(name)

        if name in self._samples:
            raise ValueError('Sample of name "%s" already exists' % name)

        self._samples[name] = sample

        if lattice is not None:
            self._set_lattice(sample, lattice)

        if select:
            self._sample = sample
            self._re_init()

        return sample

    def _set_lattice(self, sample, lattice):
        if not isinstance(lattice, HklModule.Lattice):
            a, b, c = lattice[:3]
            angles = lattice[3:]
            if self._degrees:
                angles = [np.radians(angle) for angle in angles]

            lattice = HklModule.Lattice.new(a, b, c, angles[0], angles[1], angles[2],
                                            UserUnits)

        sample.lattice_set(lattice)

    @property
    def lattice(self):
        lattice = self._sample.lattice_get(UserUnits)
        a, b, c = lattice[:3]
        angles = lattice[3:]
        if self._degrees:
            angles = [np.degrees(angle) for angle in angles]

        return a, b, c, angles[0], angles[1], angles[2]

    @lattice.setter
    def lattice(self, lattice):
        self._set_lattice(self._sample, lattice)

    @property
    def U(self):
        '''
        The U matrix
        '''
        return self._sample.U_get()

    @U.setter
    def U(self, new_u):
        self._sample.U_set(new_u)

    @property
    def ux(self):
        return self._sample.ux_get()

    @property
    def uy(self):
        return self._sample.uy_get()

    @property
    def uz(self):
        return self._sample.uz_get()

    @property
    def UB(self):
        '''
        The UB matrix
        '''
        return self._sample.UB_get()

    @UB.setter
    def UB(self, new_ub):
        self._sample.UB_set(new_ub)

    @property
    def reciprocal(self):
        '''
        The reciprocal lattice
        '''
        lattice = self._sample.lattice_get(UserUnits)
        reciprocal = lattice.copy()
        lattice.reciprocal(reciprocal)
        return reciprocal

    @property
    def reflections(self):
        '''
        All reflections for the current sample in the form:
            [(h, k, l), ...]
        '''
        return [refl.hkl_get() for refl in self._sample.reflections_get()]

    @reflections.setter
    def reflections(self, refls):
        self.clear_reflections()
        for refl in refls:
            self.add_reflection(*refl)

    def add_reflection(self, h, k, l, detector=None):
        '''
        Add a reflection, optionally specifying the detector to use
        '''
        if detector is None:
            detector = self._detector

        return self._sample.add_reflection(self._geometry, detector, h, k, l)

    def remove_reflection(self, refl):
        '''
        Remove a specific reflection
        '''
        if not isinstance(refl, HklModule.SampleReflection):
            index = self.reflections.index(refl)
            refl = self._sample.reflections_get()[index]

        return self._sample.del_reflection(refl)

    def clear_reflections(self):
        '''
        Clear all reflections for the current sample
        '''
        reflections = self._sample.reflections_get()
        for refl in reflections:
            self._sample.del_reflection(refl)

    def _refl_matrix(self, fcn):
        '''
        Get a reflection angle matrix
        '''
        sample = self._sample
        refl = sample.reflections_get()
        refl_matrix = np.zeros((len(refl), len(refl)))

        for i, r1 in enumerate(refl):
            for j, r2 in enumerate(refl):
                if i != j:
                    refl_matrix[i, j] = fcn(r1, r2)

        return refl_matrix

    @property
    def reflection_measured_angles(self):
        # TODO: typo bug report (mesured)
        return self._refl_matrix(self._sample.get_reflection_mesured_angle)

    @property
    def reflection_theoretical_angles(self):
        return self._refl_matrix(self._sample.get_reflection_theoretical_angle)

    def _re_init(self):
        if self._engine is None:
            return

        if self._geometry is None or self._detector is None or self._sample is None:
            # raise ValueError('Not all parameters set (geometry, detector, sample)')
            pass
        else:
            self._engine_list.init(self._geometry, self._detector, self._sample)

    @property
    def engines(self):
        return dict((engine.name_get(), engine)
                    for engine in self._engine_list.engines_get())

    @property
    def axes(self):
        return self._geometry.axes_names_get()

    def __getitem__(self, axis):
        return self._geometry.axis_get(axis)

    def get_axis_limits(self, axis):
        return self[axis].min_max_get(UserUnits)

    def set_axis_limits(self, axis, low, high):
        self[axis].min_max_set(low, high, UserUnits)

    def select_solution(self, sol):
        if self._solutions is None:
            raise RuntimeError('No calculation in progress')

        engine, user_sol, solutions = self._solutions
        if sol in user_sol:
            idx = user_sol.index(sol)
            sol = solutions[idx]

        engine.select_solution(sol)
        self._solutions = None

    def calc(self, h, k, l, engine=None,
             use_first=True):
        with self.using_engine(engine):
            if self.engine is None:
                raise ValueError('Engine unset')

            engine = self.engine
            try:
                solutions = self.engine.pseudo_axes_values_set([h, k, l],
                                                               UserUnits)
            except GLib.GError as ex:
                raise ValueError('Calculation failed (%s)' % ex)

            if use_first:
                # just use the first solution
                sol = solutions.items()[0]
                return [sol.geometry_get().axes_values_get(UserUnits)]
            else:
                ret = [sol.geometry_get().axes_values_get(UserUnits)
                       for sol in solutions.items()]

                if len(ret) > 1:
                    self._solutions = (engine, ret, solutions)

                return ret

    def using_engine(self, engine):
        return UsingEngine(self, engine)

    def __call__(self, start, end=None, n=100, engine=None,
                 **kwargs):
        start = np.array(start)

        # TODO better interpretations of input
        if end is not None:
            end = np.array(end)

            if start.size == end.size == 3:
                # start= [h1, k1, l1]
                # end  = [h2, k2, l2]

                # from start to end, in a linear path
                hs = np.linspace(start[0], end[0], n + 1)
                ks = np.linspace(start[1], end[1], n + 1)
                ls = np.linspace(start[2], end[2], n + 1)

            else:
                raise ValueError('Invalid start/end position')

        else:
            hkls = np.array(start)
            if hkls.ndim == 1 and hkls.size == 3:
                # single h, k, l position
                hs = [hkls[0]]
                ks = [hkls[1]]
                ls = [hkls[2]]
            elif (hkls.ndim == 2) and (3 in hkls.shape):
                if hkls.shape[0] == 3 or hkls.ndim == 1:
                    # [[h, k, l], [h, k, l], ...]
                    hs = hkls[:, 0]
                    ks = hkls[:, 1]
                    ls = hkls[:, 2]
                else:
                    # [[h, h, h, ...], [k, k, k, ...], [l, l, l, ...]]
                    hs = hkls[0, :]
                    ks = hkls[1, :]
                    ls = hkls[2, :]
            else:
                raise ValueError('Invalid set of h, k, l positions')

        with self.using_engine(engine):
            for h, k, l in zip(hs, ks, ls):
                print('calc with', h, k, l)
                yield self.calc(h, k, l, engine=None,
                                **kwargs)


class Diffractometer(PseudoPositioner):
    def __init__(self, hkl_calc, **kwargs):
        PseudoPositioner.__init__(self,
                                  forward=self.hkl_to_real,
                                  reverse=self.real_to_hkl,
                                  pseudo=['h', 'k', 'l'],
                                  **kwargs)

        self._hkl = hkl_calc

    def hkl_to_real(self, h=0.0, k=0.0, l=0.0):
        return [0, 0, 0]

    def real_to_hkl(self, **todo):
        pass
