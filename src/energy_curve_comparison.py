from math import sqrt

def root_mean_square(energy_curve_original, energy_curve_new):
    deviation = 1.20
    original_ep = [max(energy_curve_original), energy_curve_original[-1]]
    modified_ep = [max(energy_curve_new), energy_curve_new[-1]]
    border_curve = [max(energy_curve_original) * deviation, energy_curve_original[-1] * deviation]

    def rmsd(curve1, curve2):
        difference = []
        zip_object = zip(curve1, curve2)
        for list1_i, list2_i in zip_object:
            difference.append((list1_i - list2_i) ** 2)
        return sqrt((sum(difference) / 2))
    actual_rmsd = rmsd(original_ep, modified_ep)
    borderline_rmsd = rmsd(original_ep, border_curve)
    # good rmsd returns 1, bad rmsd returns 0
    if actual_rmsd < borderline_rmsd:
        return 1
    else:
        return 0
