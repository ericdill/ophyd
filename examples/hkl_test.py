from ophyd.utils.hkl import RecipCalc, hkl_module


def test():
    k6c = RecipCalc('K6C', engine='hkl')

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
    print(k6c.engine)

    return k6c

if __name__ == '__main__':
    k6c = test()
