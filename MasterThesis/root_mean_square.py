from math import sqrt


coords = [0.000000, 6.646982, 26.498723, 56.922533, 96.617970, 146.207023, 121.386505, 114.406735]

coords2 = [0.000000, 6.790366, 27.059802, 58.906111, 100.867318, 149.478833, 122.747973, 114.464051]

# Make lists same length
while len(coords) > len(coords2):
    coords2.append(coords2[-1])
while len(coords) < len(coords2):
    coords.append(coords[-1])


difference = []
zip_object = zip(coords, coords2)
for list1_i, list2_i in zip_object:
    difference.append((list1_i-list2_i)**2)

val = sqrt((sum(difference)/len(coords)))
print(val)



