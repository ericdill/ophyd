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
    from gi.repository import Hkl as hkl_module
    from gi.repository import GLib
except ImportError as ex:
    print('[!!] Failed to import Hkl library; diffractometer support'
          ' disabled (%s)' % ex,
          file=sys.stderr)

    hkl_module = None

from ..controls import PseudoPositioner


def new_detector(dtype=0):
    '''
    Create a new HKL-library detector
    '''
    return hkl_module.Detector.factory_new(hkl_module.DetectorType(dtype))


if hkl_module:
    DIFF_TYPES = tuple(sorted(hkl_module.factories().keys()))
    UserUnits = hkl_module.UnitEnum.USER
    DefaultUnits = hkl_module.UnitEnum.DEFAULT

    VALID_UNITS = (UserUnits, DefaultUnits)
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


def hkl_matrix_to_numpy(m):
    if isinstance(m, np.ndarray):
        return m

    ret = np.zeros((3, 3))
    for i in range(3):
        for j in range(3):
            ret[i, j] = m.get(i, j)

    return ret


def numpy_to_hkl_matrix(m):
    if isinstance(m, hkl_module.Matrix):
        return m

    m = np.array(m)

    hklm = hkl_euler_matrix(0, 0, 0)
    hklm.init(*m.flatten())
    return hklm


def hkl_euler_matrix(euler_x, euler_y, euler_z):
    return hkl_module.Matrix.new_euler(euler_x, euler_y, euler_z)


class HklSample(object):
    def __init__(self, calc, sample=None, units=UserUnits,
                 **kwargs):
        if sample is None:
            sample = hkl_module.Sample.new('')

        self._calc = calc
        self._sample = sample
        self._sample_dict = calc._samples
        self._units = units

        assert units in VALID_UNITS

        for name in ('lattice', 'name', 'U', 'UB', 'ux', 'uy', 'uz',
                     'reflections', ):
            value = kwargs.pop(name, None)
            if value is not None:
                try:
                    setattr(self, name, value)
                except Exception as ex:
                    ex.msg = '%s (attribute=%s)' % (ex, name)
                    raise ex(ex.msg)

        if kwargs:
            raise ValueError('Unsupported kwargs for HklSample: %s' %
                             tuple(kwargs.keys()))

    @property
    def hkl_calc(self):
        '''
        The HklCalc instance associated with the sample
        '''
        return self._calc

    @property
    def hkl_sample(self):
        '''
        The HKL library sample object
        '''
        return self._sample

    @property
    def name(self):
        '''
        The name of the currently selected sample
        '''
        return self._sample.name_get()

    @name.setter
    def name(self, new_name):
        sample = self._sample
        current = sample.name_get()
        if new_name in self._sample_dict:
            raise ValueError('Sample with that name already exists')

        sample.name_set(new_name)

        del self._sample_dict[current]
        self._sample_dict[new_name] = self

    def _set_lattice(self, sample, lattice):
        if not isinstance(lattice, hkl_module.Lattice):
            a, b, c, alpha, beta, gamma = lattice

            lattice = hkl_module.Lattice.new(a, b, c, alpha, beta, gamma,
                                             self._units)

        sample.lattice_set(lattice)

    @property
    def reciprocal(self):
        '''
        The reciprocal lattice
        '''
        lattice = self._sample.lattice_get()
        reciprocal = lattice.copy()
        lattice.reciprocal(reciprocal)
        return reciprocal.get(self._units)

    @property
    def lattice(self):
        '''
        The lattice
        '''
        lattice = self._sample.lattice_get()
        lattice = lattice.get(self._units)

        a, b, c, alpha, beta, gamma = lattice
        return a, b, c, alpha, beta, gamma

    @lattice.setter
    def lattice(self, lattice):
        self._set_lattice(self._sample, lattice)

    @property
    def U(self):
        '''
        The U matrix
        '''
        return hkl_matrix_to_numpy(self._sample.U_get())

    @U.setter
    def U(self, new_u):
        self._sample.U_set(numpy_to_hkl_matrix(new_u))

    def _get_hkl_parameter(self, param):
        return HklParameter(param, units=self._units)

    @property
    def ux(self):
        return self._get_hkl_parameter(self._sample.ux_get())

    @property
    def uy(self):
        return self._get_hkl_parameter(self._sample.uy_get())

    @property
    def uz(self):
        return self._get_hkl_parameter(self._sample.uz_get())

    @property
    def UB(self):
        '''
        The UB matrix
        '''
        return hkl_matrix_to_numpy(self._sample.UB_get())

    @UB.setter
    def UB(self, new_ub):
        self._sample.UB_set(numpy_to_hkl_matrix(new_ub))

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
            detector = self._calc._detector

        return self._sample.add_reflection(self._calc._geometry, detector,
                                           h, k, l)

    def remove_reflection(self, refl):
        '''
        Remove a specific reflection
        '''
        if not isinstance(refl, hkl_module.SampleReflection):
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

    @property
    def affine(self):
        return self._sample.affine()

    def _repr_info(self):
        repr = ['name={!r}'.format(self.name),
                'lattice={!r}'.format(self.lattice),
                'ux={!r}'.format(self.ux),
                'uy={!r}'.format(self.uy),
                'uz={!r}'.format(self.uz),
                'U={!r}'.format(self.U),
                'UB={!r}'.format(self.UB),
                'reflections={!r}'.format(self.reflections),
                ]

        return repr

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join(self._repr_info()))

    def __str__(self):
        info = self._repr_info()
        info.append('reflection_measured_angles={!r}'.format(self.reflection_measured_angles))
        info.append('reflection_theoretical_angles={!r}'.format(self.reflection_theoretical_angles))
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join(info))


class HklParameter(object):
    def __init__(self, param, units=UserUnits):
        self._param = param
        self._units = units

    @property
    def hkl_parameter(self):
        '''
        The HKL library parameter object
        '''
        return self._param

    @property
    def units(self):
        return self._units

    @units.setter
    def units(self, units):
        assert units in VALID_UNITS
        self._units = units

    @property
    def name(self):
        return self._param.name_get()

    @property
    def value(self):
        return self._param.value_get(self._units)

    @property
    def user_units(self):
        '''
        A string representing the user unit type
        '''
        return self._param.user_unit_get()

    @property
    def default_units(self):
        '''
        A string representing the default unit type
        '''
        return self._param.default_unit_get()

    @value.setter
    def value(self, value):
        self._param.value_set(value, self._units)

    @property
    def fit(self):
        return self._param.fit_get()

    @fit.setter
    def fit(self, fit):
        self._param.fit_set(fit)

    @property
    def limits(self):
        return self._param.min_max_get(self._units)

    @limits.setter
    def limits(self, (low, high)):
        self._param.min_max_set(low, high, self._units)

    def _repr_info(self):
        repr = ['name={!r}'.format(self.name),
                'limits={!r}'.format(self.limits),
                'value={!r}'.format(self.value),
                'fit={!r}'.format(self.fit),
                ]

        if self._units is UserUnits:
            repr.append('units={!r}'.format(self.user_units))
        else:
            repr.append('units={!r}'.format(self.default_units))

        return repr

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join(self._repr_info()))

    def __str__(self):
        info = self._repr_info()
        # info.append(self.)
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join(info))


class HklCalc(object):
    def __init__(self, dtype, engine='hkl',
                 sample='main', lattice=None,
                 degrees=True, units=UserUnits):
        self._engine = None  # set below with property
        self._detector = new_detector()
        self._degrees = bool(degrees)
        self._sample = None
        self._samples = {}
        self._units = units

        try:
            self._factory = hkl_module.factories()[dtype]
        except KeyError:
            raise ValueError('Invalid diffractometer type %r; '
                             'choose from: %s' % (dtype, ', '.join(DIFF_TYPES)))

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

        if isinstance(engine, hkl_module.Engine):
            self._engine = engine
        else:
            engines = self.engines
            try:
                self._engine = engines[engine]
            except KeyError:
                raise ValueError('Unknown engine name or type')

        self._re_init()

    def _get_sample(self, name):
        if isinstance(name, hkl_module.Sample):
            return name

        return self._samples[name]

    @property
    def sample_name(self):
        '''
        The name of the currently selected sample
        '''
        return self._sample.name_get()

    @sample_name.setter
    def sample_name(self, new_name):
        sample = self._sample
        sample.name = new_name

    @property
    def sample(self):
        return self._sample

    @sample.setter
    def sample(self, sample):
        if sample is self._sample:
            return
        elif sample == self._sample.name:
            return

        if isinstance(sample, HklSample):
            if sample not in self._samples.values():
                self.add_sample(sample, select=False)
        elif sample in self._samples:
            name = sample
            sample = self._samples[name]
        else:
            raise ValueError('Unknown sample type (expected HklSample)')

        self._sample = sample
        self._re_init()

    def add_sample(self, name, select=True,
                   **kwargs):
        if isinstance(name, hkl_module.Sample):
            sample = HklSample(self, name, units=self._units,
                               **kwargs)
        elif isinstance(name, HklSample):
            sample = name
        else:
            sample = HklSample(self, sample=hkl_module.Sample.new(name),
                               units=self._units,
                               **kwargs)

        if sample.name in self._samples:
            raise ValueError('Sample of name "%s" already exists' % name)

        self._samples[sample.name] = sample
        if select:
            self._sample = sample
            self._re_init()

        return sample

    def _re_init(self):
        if self._engine is None:
            return

        if self._geometry is None or self._detector is None or self._sample is None:
            # raise ValueError('Not all parameters set (geometry, detector, sample)')
            pass
        else:
            self._engine_list.init(self._geometry, self._detector,
                                   self._sample.hkl_sample)

    @property
    def engines(self):
        return dict((engine.name_get(), engine)
                    for engine in self._engine_list.engines_get())

    @property
    def axis_names(self):
        return self._geometry.axes_names_get()

    def _get_hkl_parameter(self, param):
        return HklParameter(param, units=self._units)

    def __getitem__(self, axis):
        return self._get_hkl_parameter(self._geometry.axis_get(axis))

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
                                                               self._units)
            except GLib.GError as ex:
                raise ValueError('Calculation failed (%s)' % ex)

            if use_first:
                # just use the first solution
                sol = solutions.items()[0]
                return [sol.geometry_get().axes_values_get(self._units)]
            else:
                ret = [sol.geometry_get().axes_values_get(self._units)
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
