from __future__ import print_function
from ophyd.utils.hkl import (CalcRecip, CalcE4CH, CalcK6C,
                             DiffE4CH)
import ophyd.utils.hkl as hkl_module
from ophyd.controls.positioner import Positioner
from pprint import pprint


class DumbPositioner(Positioner):
    def move(self, position, **kwargs):
        self._set_position(position)

        self._started_moving = True
        self._done_moving()

        Positioner.move(self, position, **kwargs)

    @property
    def moving(self):
        return False


def test():
    print('Diffractometer types: %s' % ', '.join(hkl_module.DIFF_TYPES))
    k6c = CalcK6C(engine='hkl')
    # or equivalently:
    # k6c = CalcRecip('K6C', engine='hkl')

    print(k6c.engines)
    print(k6c['mu'])
    print(k6c[k6c.physical_axis_names[0]].limits)
    # geometry holds physical motor information
    print('physical axes (depends on diffr. type)', k6c.physical_axis_names)
    # engine holds pseudo motor information
    print('pseudo axes (depends on engine)', k6c.pseudo_axis_names)
    print('engine parameters', k6c.parameters)
    print('1, 1, 1 -> ', list(k6c([1, 0.99, 1])))

    sample = k6c.sample
    refl = sample.add_reflection(1, 1, 1)
    sample.remove_reflection(refl)
    sample.clear_reflections()

    lim = (0.0, 20.0)
    k6c['mu'].limits = lim
    print('mu limits', k6c['mu'].limits)
    assert(k6c['mu'].limits == lim)

    k6c['h'] = 1.0
    k6c['mu'] = 0.55
    print('pseudo=', dict(k6c.engine.pseudo_axes), 'physical motors=', dict(k6c.physical_axes))

    sample.add_reflection(1, 1, 1)
    sample.add_reflection(1, 0, 1)
    sample.add_reflection(1, 0, 0)
    print(sample.reflection_measured_angles)
    print(sample.reflection_theoretical_angles)
    print(sample.reflections)

    k6c.sample.name = 'main_sample'

    sample2 = k6c.add_sample('sample2')
    try:
        k6c.add_sample('sample2')
    except ValueError:
        pass
    else:
        sample2
        raise Exception

    k6c.sample = 'main_sample'

    sample.U = [[1, 1, 1], [1, 0, 0], [1, 1, 0]]
    print('U=%s' % sample.U)
    # sample.UB = [[1, 1, 1], [1, 0, 0], [1, 1, 0]]
    print('UB=%s' % sample.UB)
    print('ux, uy, uz=%s, %s, %s' % (sample.ux, sample.uy, sample.uz))
    print('lattice=%s reciprocal=%s' % (sample.lattice, sample.reciprocal))
    print('main_sample=%s' % sample)
    # print(k6c)
    print()
    print('current engine is', k6c.engine)
    print('available engines', end=': ')
    pprint(k6c.engines)

    # TODO compute_UB affects sample state?
    # print('computed ub=%s' % sample.compute_UB([1, 1, 1], [1, 0, 1]))

    print('wavelength is %s nm (energy=%s keV)' % (k6c.wavelength, k6c.energy))

    print('hkl mode is %s (can be: %s)' % (k6c.engine.mode, k6c.engine.modes))
    print('* single position')
    list(k6c([0, 1, 0]))

    print('* 10 positions between two hkls')
    for solutions in k6c([0, 1, 0], [0, 1, 0.1], n=10):
        print('choosing', solutions[0], 'of %d solutions' % len(solutions))
        solutions[0].select()

    print('* 3 specific hkls')
    list(k6c([[0, 1, 0], [0, 1, 0.01], [0, 1, 0.02]]))

    q2_recip = CalcRecip('K6C', engine='q2')
    print('q is', q2_recip['q'])
    print('alpha is', q2_recip['alpha'])
    assert(len(list(q2_recip([[1, 2], ]))) == 1)
    assert(len(list(q2_recip([[1, 2], [3, 4]]))) == 2)
    assert(len(list(q2_recip([1, 2], [3, 4], n=20))) == 21)

    e4ch = CalcE4CH()
    print('e4ch axes:', e4ch.pseudo_axis_names, e4ch.physical_axis_names)

    positioners = [DumbPositioner(name='%s' % name) for name in
                   e4ch.physical_axis_names]

    for i, pos in enumerate(positioners):
        pos._position = 0.1 * (i + 1)

    diffr = DiffE4CH(positioners, name='my_diffractometer',
                     energy=8.0,
                     )

    calc = diffr.calc
    sample = calc.sample
    sample.add_reflection(1, 1, 1)

    pos0 = positioners[0]
    # this will run the callbacks to force a readback pseudo position calculation:
    # (not normally used, since they should be tied to real motors)
    pos0._set_position(pos0.position)

    def show_pos():
        _pseudos = [(pos.name, pos.position) for pos in diffr.pseudos.values()]
        _reals = [(pos.name, pos.position) for pos in diffr.reals.values()]

        print('pseudo positioner is at %s' % (_pseudos, ))
        print('real positioners %s' % (_reals, ))

    show_pos()
    print()
    diffr.move((1, 0, 1, 8.0), wait=True)
    print()
    show_pos()
    return k6c, diffr

if __name__ == '__main__':
    k6c, diffr = test()
