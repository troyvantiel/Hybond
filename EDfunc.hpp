

//definitions of the methods for electron density go here






//Constants for the promolecular approach of calculating the electron density and the Reduced Energy Gradient

static const int maxAtom = 18; // this code can only go up to argon
const double c1[maxAtom] = {0.2815, 2.437, 11.84, 31.34, 67.82, 120.2, 190.9,289.5,  406.3, 561.3, 760.8, 1016, 1319, 1658,2042, 2501, 3024, 3625};
const double c2[maxAtom] = {0,0, 0.06332, 0.3694, 0.8527, 1.172, 2.247,2.879, 3.049,6.984,22.42,37.17,57.95, 87.16,115.7, 158.0,   205.5,  260.0};
const double c3[maxAtom] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0,0.06358, 0.3331, 0.8878, 0.7888, 1.465, 2.170,3.369, 5.211};
const double zeta1[maxAtom] = {0.5288, 0.3379, 0.1912, 0.1390, 0.1059, 0.0884,0.0767, 0.0669, 0.0608, 0.0549, 0.0496, 0.0449,0.0411, 0.0382, 0.0358, 0.0335, 0.0315, 0.0296};
const double zeta2[maxAtom] = {1, 1, 0.9992, 0.6945, 0.5300, 0.5480,0.4532, 0.3974, 0.3994, 0.3447, 0.2511, 0.2150,0.1874, 0.1654, 0.1509, 0.1369, 0.1259, 0.1168};
const double zeta3[maxAtom] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1.0236, 0.7753, 0.5962, 0.6995, 0.5851, 0.5149,0.4974, 0.4412};
