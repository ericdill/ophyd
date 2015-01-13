from ophyd.utils.hkl import HklCalc, hkl_module


def test():
    k6c = HklCalc('K6C')

    print(k6c.axis_names, k6c.engines)
    print(k6c['mu'])
    print(k6c[k6c.axis_names[0]].limits)
    print('1, 1, 1 -> ', list(k6c([1, 1, 1])))

    sample = k6c.sample
    refl = sample.add_reflection(1, 1, 1)
    sample.remove_reflection(refl)
    sample.clear_reflections()

    lim = (0.0, 20.0)
    k6c['mu'].limits = lim
    print('mu limits', k6c['mu'].limits)
    assert(k6c['mu'].limits == lim)

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
    return k6c

if __name__ == '__main__':
    k6c = test()
