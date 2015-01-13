from ophyd.utils.hkl import HklCalc, hkl_module


def test():
    k6c = HklCalc('K6C')

    print(k6c.axis_names, k6c.engines)
    print(k6c['mu'])
    print(k6c[k6c.axis_names[0]].limits)
    print('1, 1, 1 -> ', list(k6c([1, 1, 1])))
    refl = k6c.add_reflection(1, 1, 1)
    k6c.remove_reflection(refl)
    k6c.clear_reflections()

    lim = (0.0, 20.0)
    k6c['mu'].limits = lim
    print('mu limits', k6c['mu'].limits)
    assert(k6c['mu'].limits == lim)

    k6c.add_reflection(1, 1, 1)
    k6c.add_reflection(1, 0, 1)
    k6c.add_reflection(1, 0, 0)
    print(k6c.reflection_measured_angles)
    print(k6c.reflection_theoretical_angles)
    print(k6c.reflections)

    k6c.add_sample('sample2')
    try:
        k6c.add_sample('sample2')
    except ValueError:
        pass
    else:
        raise Exception

    k6c.sample = 'main'

    print('U=%s' % k6c.U)
    print('UB=%s' % k6c.UB)
    print('ux, uy, uz=%s, %s, %s' % (k6c.ux, k6c.uy, k6c.uz))
    return k6c

if __name__ == '__main__':
    k6c = test()
