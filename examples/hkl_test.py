from ophyd.utils.hkl import HklCalc, HklModule


def test():
    k6c = HklCalc('K6C')

    print(k6c.axes, k6c.engines)
    print(k6c['mu'])
    print(k6c.get_axis_limits(k6c.axes[0]))
    print('1, 1, 1 -> ', list(k6c([1, 1, 1])))
    refl = k6c.add_reflection(1, 1, 1)
    k6c.remove_reflection(refl)
    k6c.clear_reflections()

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

    return k6c

if __name__ == '__main__':
    k6c = test()
    print('hkl is', HklModule)
