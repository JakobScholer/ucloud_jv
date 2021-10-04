from math import sqrt


def root_mean_square(energy_curve_original, energy_curve_new):
    c1 = [max(energy_curve_original), energy_curve_original[-1]]
    c2 = [max(energy_curve_new), energy_curve_new[-1]]

    difference = []
    zip_object = zip(c1, c2)
    for list1_i, list2_i in zip_object:
        difference.append((list1_i-list2_i)**2)
    return sqrt((sum(difference)/2))




