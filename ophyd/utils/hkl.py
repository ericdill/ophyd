# vi: ts=4 sw=4 sts=4 expandtab
'''
:mod:`ophyd.utils.hkl` - HKL calculation utilities
==================================================

.. module:: ophyd.utils.hkl
   :synopsis:

'''

from __future__ import print_function
import inspect
import sys
from collections import OrderedDict

import numpy as np

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
                    ex.message = '%s (attribute=%s)' % (ex, name)
                    raise ex(ex.message)

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

        # TODO: notes mention that lattice should not change, but is it alright
        #       if init() is called again? or should reflections be cleared,
        #       etc?

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

    def _get_parameter(self, param):
        return Parameter(param, units=self._units)

    @property
    def ux(self):
        return self._get_parameter(self._sample.ux_get())

    @property
    def uy(self):
        return self._get_parameter(self._sample.uy_get())

    @property
    def uz(self):
        return self._get_parameter(self._sample.uz_get())

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


class Parameter(object):
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


class Solution(object):
    def __init__(self, engine, list_item):
        self._list_item = list_item.copy()
        self._geometry = list_item.geometry_get().copy()
        self._engine = engine

    def __getitem__(self, axis):
        return self._geometry.axis_get(axis)

    @property
    def axis_names(self):
        return self._geometry.axes_names_get()

    @property
    def axis_values(self):
        return self._geometry.axes_values_get(self.units)

    @property
    def units(self):
        return self._engine.units

    def set_wavelength(self, wavelength):
        # TODO
        self._geometry.wavelength_set(wavelength)

    def select(self):
        self._engine._engine_list.select_solution(self._list_item)

    def _repr_info(self):
        repr = ['{!r}'.format(self.axis_values),
                'units={!r}'.format(self.units),
                ]

        return repr

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join(self._repr_info()))


class Engine(object):
    def __init__(self, calc, engine, engine_list):
        self._calc = calc
        self._engine = engine
        self._engine_list = engine_list
        self._solutions = None

    @property
    def name(self):
        return self._engine.name_get()

    @property
    def mode(self):
        '''
        HKL calculation mode (see also `HklCalc.modes`)
        '''
        return self._engine.current_mode_get()

    @mode.setter
    def mode(self, mode):
        if mode not in self.modes:
            raise ValueError('Unrecognized mode %r; '
                             'choose from: %s' % (mode, ', '.join(self.modes))
                             )

        return self._engine.current_mode_set(mode)

    @property
    def modes(self):
        return self._engine.modes_names_get()

    @property
    def solutions(self):
        return tuple(self._solutions)

    @property
    def parameters(self):
        return self._engine.parameters_names_get()

    @property
    def pseudo_axis_names(self):
        return self._engine.pseudo_axes_names_get()

    @property
    def pseudo_axis_values(self):
        return self._engine.pseudo_axes_values_get(self.units)

    @property
    def pseudo_axes(self):
        keys = self.pseudo_axis_names
        values = self.pseudo_axis_values
        return OrderedDict(zip(keys, values))

    @pseudo_axis_values.setter
    def pseudo_axis_values(self, values):
        try:
            geometry_list = self._engine.pseudo_axes_values_set(values, self.units)
        except GLib.GError as ex:
            raise ValueError('Calculation failed (%s)' % ex)

        self._solutions = [Solution(self, item)
                           for item in geometry_list.items()]

    def __getitem__(self, name):
        try:
            return self.pseudo_axes[name]
        except KeyError:
            raise ValueError('Unknown axis name: %s' % name)

    def __setitem__(self, name, value):
        values = self.pseudo_axis_values
        try:
            idx = self.pseudo_axis_names.index(name)
        except IndexError:
            raise ValueError('Unknown axis name: %s' % name)

        values[idx] = float(value)
        self.pseudo_axis_values = values

    @property
    def units(self):
        return self._calc._units

    @property
    def engine(self):
        return self._engine

    def _repr_info(self):
        repr = ['parameters={!r}'.format(self.parameters),
                'pseudo_axes={!r}'.format(dict(self.pseudo_axes)),
                'mode={!r}'.format(self.mode),
                'modes={!r}'.format(self.modes),
                'units={!r}'.format(self.units),
                ]

        return repr

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join(self._repr_info()))


class CalcRecip(object):
    def __init__(self, dtype, engine='hkl',
                 sample='main', lattice=None,
                 degrees=True, units=UserUnits,
                 lock_engine=False):
        self._engine = None  # set below with property
        self._detector = new_detector()
        self._degrees = bool(degrees)
        self._sample = None
        self._samples = {}
        self._units = units
        self._lock_engine = bool(lock_engine)

        try:
            self._factory = hkl_module.factories()[dtype]
        except KeyError:
            raise ValueError('Invalid diffractometer type %r; '
                             'choose from: %s' % (dtype, ', '.join(DIFF_TYPES)))

        self._geometry = self._factory.create_new_geometry()
        self._engine_list = self._factory.create_new_engine_list()

        if sample is not None:
            self.add_sample(sample, lattice=lattice)

        self.engine = engine

    @property
    def engine_locked(self):
        '''
        If set, do not allow the engine to be changed post-initialization
        '''
        return self._lock_engine

    @property
    def engine(self):
        return self._engine

    @engine.setter
    def engine(self, engine):
        if engine is self._engine:
            return

        if self._lock_engine and self._engine is not None:
            raise ValueError('Engine is locked on this %s instance' %
                             self.__class__.__name__)

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
        return dict((engine.name_get(), Engine(self, engine, self._engine_list))
                    for engine in self._engine_list.engines_get())

    @property
    def parameters(self):
        return self._engine.parameters

    @property
    def physical_axis_names(self):
        return self._geometry.axes_names_get()

    @property
    def physical_axis_values(self):
        return [self[name].value
                for name in self._geometry.axes_names_get()]

    @property
    def physical_axes(self):
        keys = self.physical_axis_names
        values = self.physical_axis_values
        return OrderedDict(zip(keys, values))

    @property
    def pseudo_axis_names(self):
        '''Pseudo axis names from the current engine'''
        return self._engine.pseudo_axis_names

    @property
    def pseudo_axis_values(self):
        '''Pseudo axis positions/values from the current engine'''
        return self._engine.pseudo_axis_values

    @property
    def pseudo_axes(self):
        '''Dictionary of axis name to position'''
        return self._engine.pseudo_axes

    def _get_parameter(self, param):
        return Parameter(param, units=self._units)

    def __getitem__(self, axis):
        if axis in self.physical_axis_names:
            return self._get_parameter(self._geometry.axis_get(axis))
        elif axis in self.pseudo_axis_names:
            return self._engine[axis]

    def __setitem__(self, axis, value):
        if axis in self.physical_axis_names:
            param = self[axis]
            param.value = value
        elif axis in self.pseudo_axis_names:
            self._engine[axis] = value

    def calc(self, position, engine=None,
             use_first=False):
        # TODO default should probably not be `use_first` (or remove
        # completely?)
        with self.using_engine(engine):
            if self.engine is None:
                raise ValueError('Engine unset')

            engine = self.engine
            self.engine.pseudo_axis_values = position

            solutions = self.engine.solutions

            if use_first:
                # just use the first solution
                solutions[0].select()

            return solutions

    def using_engine(self, engine):
        return UsingEngine(self, engine)

    def calc_linear_path(self, start, end, n, num_params=0, **kwargs):
        # start = [h1, k1, l1]
        # end   = [h2, k2, l2]

        # from start to end, in a linear path
        singles = [np.linspace(start[i], end[i], n + 1)
                   for i in range(num_params)]

        return list(zip(*singles))

    def _get_path_fcn(self, path_type):
        try:
            return getattr(self, 'calc_%s_path' % (path_type))
        except AttributeError:
            raise ValueError('Invalid path type specified (%s)' % path_type)

    def get_path(self, start, end=None, n=100,
                 path_type='linear', **kwargs):
        num_params = len(self.pseudo_axis_names)

        start = np.array(start)

        path_fcn = self._get_path_fcn(path_type)

        if end is not None:
            end = np.array(end)
            if start.size == end.size == num_params:
                return path_fcn(start, end, n, num_params=num_params,
                                **kwargs)

        else:
            positions = np.array(start)
            if positions.ndim == 1 and positions.size == num_params:
                # single position
                return [list(positions)]
            elif positions.ndim == 2:
                if positions.shape[0] == 1 and positions.size == num_params:
                    # [[h, k, l], ]
                    return [positions[0]]
                elif positions.shape[0] == num_params:
                    # [[h, k, l], [h, k, l], ...]
                    return [positions[i, :] for i in range(num_params)]

        raise ValueError('Invalid set of %s positions' %
                         ', '.join(self.pseudo_axis_names))

    def __call__(self, start, end=None, n=100, engine=None,
                 path_type='linear', **kwargs):

        with self.using_engine(engine):
            for pos in self.get_path(start, end=end, n=n,
                                     path_type=path_type, **kwargs):
                yield self.calc(pos, engine=None,
                                **kwargs)

    def _repr_info(self):
        repr = ['engine={!r}'.format(self.engine.name),
                'detector={!r}'.format(self._detector),
                'sample={!r}'.format(self._sample),
                'samples={!r}'.format(self._samples),
                ]

        return repr

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join(self._repr_info()))

    def __str__(self):
        info = self._repr_info()
        return '{}({})'.format(self.__class__.__name__,
                               ', '.join(info))


class Diffractometer(PseudoPositioner):
    def __init__(self, calc_class, calc_kw=None,
                 **kwargs):

        if isinstance(calc_class, CalcRecip):
            self._calc = calc_class
        elif inspect.isclass(calc_class):
            if calc_kw is None:
                calc_kw = {}

            calc_kw = dict(calc_kw)
            calc_kw.pop('lock_engine', True)
            self._calc = calc_class(lock_engine=True, **calc_kw)
        else:
            raise ValueError('Must specify a calculation class or an instance '
                             'of CalcRecip (or a derived class). Got: %s' %
                             calc_class)

        if not self._calc.engine_locked:
            # Reason for this is that the engine determines the pseudomotor
            # names, so if the engine is switched from underneath, the
            # pseudomotor will no longer function properly
            raise ValueError('Calculation engine must have lock_engine set')

        PseudoPositioner.__init__(self,
                                  [],
                                  forward=self.hkl_to_real,
                                  reverse=self.real_to_hkl,
                                  pseudo=['h', 'k', 'l'],
                                  **kwargs)

    def hkl_to_real(self, h=0.0, k=0.0, l=0.0):
        return [0, 0, 0]

    def real_to_hkl(self, **todo):
        pass


def _create_classes(class_suffix, dtype):
    '''
    Create reciprocal calculation and diffractometer classes
    for a specific type of diffractometer.
    '''
    # - calculation
    def calc_init(self, **kwargs):
        CalcRecip.__init__(self, dtype, **kwargs)

    calc_name = 'Calc%s' % class_suffix
    _dict = dict(__init__=calc_init,
                 __doc__='Reciprocal space calculation helper for %s' % dtype)
    globals()[calc_name] = type(calc_name, (CalcRecip, ), _dict)

    calc_class = globals()[calc_name]

    # - diffractometer pseudomotor
    def diffr_init(self, **kwargs):
        Diffractometer.__init__(self, calc_class, **kwargs)

    diffr_class = 'Diff%s' % class_suffix
    _dict = dict(__init__=diffr_init,
                 __doc__='%s diffractometer pseudomotor' % dtype)
    globals()[diffr_class] = type(diffr_class, (Diffractometer, ), _dict)


# TODO can we make these names a bit better? ugh
_create_classes('E4CH', 'E4CH')
_create_classes('E4CV', 'E4CV')
_create_classes('E6C', 'E6C')
_create_classes('K4CV', 'K4CV')
_create_classes('K6C', 'K6C')
_create_classes('Petra3_p09_eh2', 'PETRA3 P09 EH2')
_create_classes('SoleilMars', 'SOLEIL MARS')
_create_classes('SoleilSiriusKappa', 'SOLEIL SIRIUS KAPPA')
_create_classes('SoleilSiriusTurret', 'SOLEIL SIRIUS TURRET')
_create_classes('SoleilSixsMed1p2', 'SOLEIL SIXS MED1+2')
_create_classes('SoleilSixsMed2p2', 'SOLEIL SIXS MED2+2')
_create_classes('SoleilSixs', 'SOLEIL SIXS')
_create_classes('Med2p3', 'MED2+3')
_create_classes('TwoC', 'TwoC')
_create_classes('Zaxis', 'ZAXIS')
