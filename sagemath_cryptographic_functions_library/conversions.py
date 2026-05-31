from sage.all import *
from .helpers import interpolate_truth_table


def univariate_to_truth_table(n, polynomial):
    r"""
    Construct the truth table of a univariate polynomial over GF(2^n).
    Defined by evaluating `f` at every element of `GF(2^n)` in order, where elements are identified with integers via their integer representation.
    
    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import univariate_to_truth_table
        sage: F.<a> = GF(2^3)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^3
        sage: tt = univariate_to_truth_table(3, polynomial); tt
        [0, 1, 3, 4, 5, 6, 7, 2]

        sage: F.<a> = GF(2^8)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = (a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^192 + (a^6 + a^3)*x^160 + (a^6 + a^4 + a^2 + 1)*x^144 + (a^5 + a^4 + a^3 + a)*x^136 + (a^6 + a^5 + a^4 + a^3 + a^2 + 1)*x^132 + (a^7 + a^6 + a^5 + a^2)*x^130 + (a^7 + a^6 + a^5 + 1)*x^129 + (a^7 + a^6 + a^4 + a^2 + a + 1)*x^128 + (a^2 + a + 1)*x^96 + (a^6 + a^3)*x^80 + (a^7 + a^6 + a^4 + a^2 + a)*x^72 + (a^7 + 1)*x^68 + (a^5 + a^4 + a + 1)*x^66 + (a^7 + a^4)*x^65 + (a^4 + a^3)*x^64 + (a^7 + a^5 + a^4 + a^3)*x^48 + (a^6 + a^5 + a^4 + a^2)*x^40 + (a^7 + a^6 + a^4 + a + 1)*x^36 + (a^5 + a^4 + 1)*x^34 + a^7*x^33 + (a^5 + a^2 + a)*x^32 + (a^7 + a^5 + a^4 + 1)*x^24 + (a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^20 + (a^5 + a^2)*x^18 + (a^7 + a^6 + a^5 + a^4 + a^3)*x^17 + (a^2 + a)*x^16 + (a^7 + a^6 + a^5 + a^4 + a^3 + 1)*x^12 + (a^5 + a^4 + a^2 + a)*x^10 + (a^6 + a^2 + a + 1)*x^9 + (a^5 + a^4 + a^3 + a^2 + 1)*x^8 + (a^5 + a^4 + a^3 + a + 1)*x^6 + (a^7 + a^5 + a^4 + 1)*x^5 + (a^7 + a^6 + a^3 + a^2 + a)*x^4 + (a^7 + a^6 + a^5 + a^2 + 1)*x^3 + (a^7 + a^4 + a^3 + a^2 + 1)*x^2 + (a^7 + a^6 + a^5 + a^4 + a^3 + a + 1)*x
        sage: tt = univariate_to_truth_table(8, polynomial); tt
        [0, 146, 38, 74, 75, 170, 36, 59, 147, 87, 114, 72, 39, 144, 143, 198, 67, 94, 207, 44, 30, 112, 219, 75, 40, 99, 99, 214, 138, 178, 136, 78, 131, 88, 158, 187, 50, 154, 102, 48, 30, 147, 196, 183, 80, 174, 195, 195, 35, 119, 148, 62, 132, 163, 122, 163, 70, 68, 54, 202, 30, 111, 39, 168, 73, 221, 157, 247, 226, 5, 127, 102, 36, 230, 55, 11, 112, 193, 42, 101, 179, 168, 205, 40, 14, 102, 57, 175, 38, 107, 159, 44, 100, 90, 148, 84, 4, 217, 235, 200, 85, 251, 243, 163, 103, 236, 79, 58, 201, 49, 168, 174, 29, 79, 88, 244, 90, 123, 86, 137, 134, 130, 4, 254, 62, 73, 245, 124, 145, 239, 189, 61, 196, 201, 161, 82, 226, 202, 9, 223, 72, 19, 234, 79, 29, 236, 155, 148, 94, 220, 145, 237, 150, 49, 215, 142, 42, 254, 34, 8, 103, 80, 112, 185, 200, 140, 150, 44, 26, 123, 202, 85, 74, 88, 211, 63, 8, 176, 181, 243, 177, 122, 69, 112, 141, 99, 247, 231, 203, 86, 248, 155, 37, 93, 251, 125, 144, 155, 7, 242, 168, 134, 177, 97, 226, 191, 178, 17, 16, 231, 100, 109, 179, 55, 142, 244, 101, 196, 214, 137, 57, 235, 195, 239, 29, 44, 248, 55, 82, 16, 254, 66, 158, 249, 188, 37, 46, 58, 69, 175, 203, 117, 132, 196, 146, 95, 148, 167, 176, 88, 56, 46, 22, 141, 215, 178]

        sage: F.<a> = GF(2^10)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = (a^7 + a^5 + a)*x^768 + (a^8 + a^6 + a^4 + a^3 + 1)*x^640 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a)*x^576 + (a^7 + a^6 + a^4 + a^2 + a + 1)*x^544 + (a^9 + a^8 + a^6 + a^5 + a^4 + a^2 + 1)*x^528 + (a^9 + a^8 + a^7 + a^5 + a^4 + a^3 + a^2 + 1)*x^520 + (a^8 + a^6 + a^3 + a)*x^516 + (a^6 + a^5 + a^4 + a^2 + a)*x^514 + (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^513 + (a^9 + a^7 + a^6 + a^3 + a^2 + a)*x^512 + (a^5 + a^4 + 1)*x^384 + (a^9 + a^7 + a^5 + a^4 + a + 1)*x^320 + (a^9 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^288 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a)*x^272 + (a^5 + a^4 + a)*x^264 + (a^8 + a^5 + a^3 + 1)*x^260 + (a^9 + a^8 + a^7 + a^4 + a^3 + 1)*x^258 + (a^7 + a^6 + a^5 + a^4 + a^2 + a + 1)*x^257 + (a^5 + a^4 + 1)*x^256 + (a^9 + a^8 + a^4 + 1)*x^192 + (a^9 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^160 + (a^8 + a^6 + a^4 + a^2 + 1)*x^144 + (a^9 + a^7 + a^6)*x^136 + (a^8 + a^7 + a^6 + a^3 + a + 1)*x^132 + (a^9 + a^8 + a^7 + a^5 + a^4 + a^3 + a^2)*x^130 + (a^6 + a^5 + a^4 + a^3 + a)*x^129 + (a^8 + a^5 + a^4 + a + 1)*x^128 + (a^9 + a^8 + a^6 + a^5 + 1)*x^96 + (a^8 + a^7 + a^5 + a + 1)*x^80 + (a^9 + a^8 + a^4 + a^3 + 1)*x^72 + (a^6 + a^5 + a^4 + a^3 + 1)*x^68 + (a^9 + a^8 + a^7 + a^5 + a^4 + a^2 + 1)*x^66 + (a^9 + a^8 + a^6 + a^5 + a^4)*x^65 + (a^9 + a^7 + a^6 + 1)*x^64 + (a^8 + a^7 + a^6 + a^5)*x^48 + (a^9 + a^7 + a^6 + a^3)*x^40 + (a^9 + a^6 + a^5 + a)*x^36 + (a^8 + a^6 + 1)*x^34 + (a^8 + a^7 + a^6 + a^5 + a^3 + a)*x^33 + (a^6 + a^4 + a^3 + 1)*x^32 + (a^7 + a^6 + a^4 + a^2 + a + 1)*x^24 + (a^8 + a^6 + a^5 + a^3 + 1)*x^20 + (a^7 + a^6 + a^4 + a^3)*x^18 + (a^9 + a^2 + a)*x^17 + (a^3 + a + 1)*x^16 + (a^8 + a^6 + a^4 + a^3 + a^2 + 1)*x^12 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a)*x^10 + (a^8 + a^5 + a^4 + a^3 + a^2 + 1)*x^9 + (a^8 + a^4 + a^3 + a^2 + a)*x^8 + (a^7 + a^4)*x^6 + (a^7 + a^6 + a^5 + a^2 + 1)*x^5 + (a^9 + a^7 + a^6 + a^4 + a^2 + 1)*x^4 + (a^5 + a^4 + a^2)*x^3 + (a^7 + a^6 + a^5 + a^3 + a)*x^2 + (a^8 + 1)*x    
        sage: tt = univariate_to_truth_table(10, polynomial); tt
        [0, 464, 928, 55, 247, 707, 110, 29, 494, 120, 561, 480, 220, 686, 58, 15, 988, 343, 240, 60, 981, 698, 960, 232, 440, 885, 747, 609, 116, 349, 30, 880, 467, 111, 596, 431, 285, 837, 419, 444, 542, 996, 486, 603, 789, 267, 980, 909, 670, 121, 405, 309, 686, 941, 668, 472, 729, 120, 429, 331, 812, 617, 865, 99, 934, 944, 222, 655, 799, 237, 862, 747, 570, 618, 317, 810, 838, 242, 888, 651, 907, 198, 127, 373, 972, 869, 769, 495, 413, 662, 534, 858, 31, 240, 173, 517, 737, 667, 446, 899, 609, 511, 519, 990, 350, 354, 638, 5, 27, 963, 2, 413, 605, 380, 398, 232, 547, 742, 713, 75, 616, 271, 452, 228, 979, 848, 838, 386, 251, 990, 215, 437, 444, 381, 681, 47, 393, 746, 474, 254, 267, 396, 609, 161, 963, 957, 867, 346, 634, 480, 483, 62, 315, 259, 484, 923, 327, 667, 673, 826, 74, 771, 65, 335, 308, 409, 518, 236, 795, 20, 879, 551, 928, 843, 237, 577, 995, 1009, 868, 305, 611, 405, 477, 108, 824, 876, 960, 467, 893, 205, 188, 331, 629, 150, 641, 549, 892, 635, 177, 1009, 885, 464, 1022, 796, 953, 760, 11, 781, 700, 772, 708, 315, 843, 279, 10, 17, 54, 456, 49, 904, 4, 542, 826, 871, 592, 223, 643, 587, 864, 523, 138, 934, 371, 954, 479, 337, 390, 171, 531, 377, 520, 988, 599, 452, 966, 502, 160, 215, 673, 819, 641, 340, 682, 220, 435, 386, 502, 124, 11, 966, 430, 960, 874, 835, 888, 692, 762, 369, 741, 205, 94, 49, 786, 451, 611, 757, 948, 641, 508, 654, 534, 129, 792, 968, 885, 518, 322, 630, 888, 670, 674, 259, 793, 283, 506, 447, 981, 629, 624, 407, 625, 53, 237, 238, 269, 944, 91, 161, 402, 203, 1021, 227, 554, 209, 771, 959, 880, 623, 352, 568, 148, 216, 433, 954, 130, 810, 670, 881, 616, 610, 818, 383, 955, 85, 472, 113, 385, 662, 40, 376, 361, 410, 1017, 333, 247, 934, 289, 55, 474, 367, 821, 455, 654, 686, 908, 491, 673, 357, 154, 281, 593, 567, 812, 269, 955, 57, 511, 58, 778, 113, 644, 952, 987, 836, 364, 948, 95, 866, 430, 212, 331, 402, 899, 285, 861, 34, 300, 20, 693, 558, 1021, 289, 335, 630, 833, 575, 354, 447, 85, 719, 349, 377, 928, 451, 75, 907, 399, 8, 197, 167, 583, 98, 22, 912, 429, 108, 177, 930, 743, 947, 352, 407, 15, 703, 128, 981, 681, 955, 148, 37, 388, 882, 544, 616, 250, 757, 783, 163, 748, 775, 411, 405, 830, 375, 369, 667, 237, 320, 791, 430, 446, 320, 689, 1004, 801, 59, 375, 904, 929, 793, 276, 15, 251, 935, 742, 772, 195, 870, 958, 440, 674, 739, 780, 680, 342, 693, 913, 465, 754, 757, 111, 698, 737, 627, 496, 193, 71, 817, 44, 703, 733, 521, 118, 257, 446, 654, 271, 129, 781, 196, 110, 516, 341, 376, 710, 782, 187, 820, 610, 78, 806, 845, 1004, 683, 248, 1016, 22, 693, 59, 223, 860, 605, 55, 881, 355, 902, 305, 403, 327, 859, 735, 644, 579, 955, 738, 349, 637, 39, 410, 391, 188, 258, 98, 923, 899, 680, 176, 988, 64, 655, 74, 194, 272, 125, 604, 374, 790, 415, 867, 941, 441, 969, 518, 561, 644, 784, 514, 465, 160, 662, 864, 785, 600, 906, 673, 308, 327, 454, 651, 77, 755, 406, 518, 804, 389, 322, 566, 182, 1012, 215, 894, 538, 29, 967, 861, 704, 855, 873, 814, 343, 853, 201, 106, 433, 474, 418, 476, 995, 444, 337, 599, 253, 561, 312, 739, 941, 861, 1014, 201, 549, 277, 602, 440, 176, 119, 961, 784, 737, 772, 854, 858, 335, 284, 748, 516, 947, 938, 958, 907, 472, 296, 666, 432, 69, 866, 820, 195, 722, 260, 752, 483, 80, 651, 667, 341, 770, 871, 910, 883, 477, 467, 734, 766, 948, 193, 110, 170, 578, 944, 251, 226, 494, 37, 1019, 154, 259, 598, 620, 464, 941, 554, 434, 746, 821, 412, 480, 613, 94, 763, 638, 712, 10, 118, 791, 892, 602, 894, 957, 818, 438, 54, 785, 835, 547, 683, 223, 747, 728, 175, 319, 982, 1, 757, 199, 714, 703, 308, 226, 562, 419, 789, 570, 985, 177, 495, 804, 538, 662, 193, 424, 114, 860, 1022, 371, 116, 190, 818, 298, 853, 778, 271, 243, 593, 490, 335, 785, 343, 334, 695, 781, 406, 619, 541, 862, 758, 498, 222, 633, 780, 1004, 1002, 751, 894, 60, 236, 525, 833, 999, 269, 16, 68, 798, 600, 161, 40, 150, 733, 902, 1003, 247, 77, 754, 578, 698, 670, 216, 859, 858, 309, 151, 969, 44, 708, 196, 894, 825, 170, 334, 553, 394, 575, 846, 849, 103, 339, 966, 772, 982, 972, 763, 733, 429, 357, 950, 845, 985, 317, 791, 223, 178, 687, 865, 116, 1021, 836, 296, 729, 754, 275, 155, 951, 120, 354, 441, 243, 623, 633, 326, 209, 425, 704, 605, 814, 500, 30, 871, 713, 1015, 256, 640, 29, 474, 741, 641, 193, 738, 296, 750, 74, 459, 776, 810, 339, 822, 708, 627, 882, 386, 486, 693, 873, 637, 837, 948, 652, 58, 418, 695, 850, 512, 567, 475, 781, 678, 491, 483, 1000, 423, 60, 918, 377, 148, 549, 619, 89, 592, 409, 625, 860, 755, 892, 880, 640, 203, 725, 379, 111, 390, 501, 447, 118, 635, 750, 605, 167, 595, 245, 930, 389, 149, 552, 733, 30, 684, 502, 743, 249, 431, 969, 77, 299, 232, 277, 373, 206, 745, 678, 356, 59, 446, 447, 409, 27, 634, 47, 240, 577, 217, 525, 310, 858, 550, 714, 595, 219, 517, 301, 592, 5, 319, 222, 752, 707, 682, 629, 959, 849, 220, 992, 392, 386, 429, 142, 258, 469, 542, 88, 301, 713, 507, 525, 156, 933, 883, 236, 479, 514, 374, 892, 427, 683, 571, 538, 88, 32, 37, 136, 302, 395, 618, 775, 259, 322, 257, 80, 432, 300, 651, 525, 788, 187, 997, 97, 668, 494, 340, 154, 453, 595, 331, 819, 392, 707, 575]
    """
    F = GF(2**n, 'a')
    return [polynomial(F.from_integer(i)).to_integer() for i in range(F.order())]


def truth_table_to_univariate(n, tt):
    r"""
    Recover the unique univariate polynomial over GF(2^n) from a truth table.
    Defined by applying Lagrange interpolation to the point set `{(i, T[i]) : i in GF(2^n)}`, where integers are identified with field elements via their integer representation.    

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``tt`` -- a truth table (look up table) represented as a list of integers

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import truth_table_to_univariate
        sage: tt = [0, 1, 3, 4, 5, 6, 7, 2]
        sage: polynomial = truth_table_to_univariate(3, tt); polynomial
        x^3

        sage: tt = [0, 146, 38, 74, 75, 170, 36, 59, 147, 87, 114, 72, 39, 144, 143, 198, 67, 94, 207, 44, 30, 112, 219, 75, 40, 99, 99, 214, 138, 178, 136, 78, 131, 88, 158, 187, 50, 154, 102, 48, 30, 147, 196, 183, 80, 174, 195, 195, 35, 119, 148, 62, 132, 163, 122, 163, 70, 68, 54, 202, 30, 111, 39, 168, 73, 221, 157, 247, 226, 5, 127, 102, 36, 230, 55, 11, 112, 193, 42, 101, 179, 168, 205, 40, 14, 102, 57, 175, 38, 107, 159, 44, 100, 90, 148, 84, 4, 217, 235, 200, 85, 251, 243, 163, 103, 236, 79, 58, 201, 49, 168, 174, 29, 79, 88, 244, 90, 123, 86, 137, 134, 130, 4, 254, 62, 73, 245, 124, 145, 239, 189, 61, 196, 201, 161, 82, 226, 202, 9, 223, 72, 19, 234, 79, 29, 236, 155, 148, 94, 220, 145, 237, 150, 49, 215, 142, 42, 254, 34, 8, 103, 80, 112, 185, 200, 140, 150, 44, 26, 123, 202, 85, 74, 88, 211, 63, 8, 176, 181, 243, 177, 122, 69, 112, 141, 99, 247, 231, 203, 86, 248, 155, 37, 93, 251, 125, 144, 155, 7, 242, 168, 134, 177, 97, 226, 191, 178, 17, 16, 231, 100, 109, 179, 55, 142, 244, 101, 196, 214, 137, 57, 235, 195, 239, 29, 44, 248, 55, 82, 16, 254, 66, 158, 249, 188, 37, 46, 58, 69, 175, 203, 117, 132, 196, 146, 95, 148, 167, 176, 88, 56, 46, 22, 141, 215, 178]
        sage: polynomial = truth_table_to_univariate(8, tt); polynomial
        (a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^192 + (a^6 + a^3)*x^160 + (a^6 + a^4 + a^2 + 1)*x^144 + (a^5 + a^4 + a^3 + a)*x^136 + (a^6 + a^5 + a^4 + a^3 + a^2 + 1)*x^132 + (a^7 + a^6 + a^5 + a^2)*x^130 + (a^7 + a^6 + a^5 + 1)*x^129 + (a^7 + a^6 + a^4 + a^2 + a + 1)*x^128 + (a^2 + a + 1)*x^96 + (a^6 + a^3)*x^80 + (a^7 + a^6 + a^4 + a^2 + a)*x^72 + (a^7 + 1)*x^68 + (a^5 + a^4 + a + 1)*x^66 + (a^7 + a^4)*x^65 + (a^4 + a^3)*x^64 + (a^7 + a^5 + a^4 + a^3)*x^48 + (a^6 + a^5 + a^4 + a^2)*x^40 + (a^7 + a^6 + a^4 + a + 1)*x^36 + (a^5 + a^4 + 1)*x^34 + a^7*x^33 + (a^5 + a^2 + a)*x^32 + (a^7 + a^5 + a^4 + 1)*x^24 + (a^7 + a^6 + a^5 + a^4 + a^2 + 1)*x^20 + (a^5 + a^2)*x^18 + (a^7 + a^6 + a^5 + a^4 + a^3)*x^17 + (a^2 + a)*x^16 + (a^7 + a^6 + a^5 + a^4 + a^3 + 1)*x^12 + (a^5 + a^4 + a^2 + a)*x^10 + (a^6 + a^2 + a + 1)*x^9 + (a^5 + a^4 + a^3 + a^2 + 1)*x^8 + (a^5 + a^4 + a^3 + a + 1)*x^6 + (a^7 + a^5 + a^4 + 1)*x^5 + (a^7 + a^6 + a^3 + a^2 + a)*x^4 + (a^7 + a^6 + a^5 + a^2 + 1)*x^3 + (a^7 + a^4 + a^3 + a^2 + 1)*x^2 + (a^7 + a^6 + a^5 + a^4 + a^3 + a + 1)*x

        sage: tt = [0, 464, 928, 55, 247, 707, 110, 29, 494, 120, 561, 480, 220, 686, 58, 15, 988, 343, 240, 60, 981, 698, 960, 232, 440, 885, 747, 609, 116, 349, 30, 880, 467, 111, 596, 431, 285, 837, 419, 444, 542, 996, 486, 603, 789, 267, 980, 909, 670, 121, 405, 309, 686, 941, 668, 472, 729, 120, 429, 331, 812, 617, 865, 99, 934, 944, 222, 655, 799, 237, 862, 747, 570, 618, 317, 810, 838, 242, 888, 651, 907, 198, 127, 373, 972, 869, 769, 495, 413, 662, 534, 858, 31, 240, 173, 517, 737, 667, 446, 899, 609, 511, 519, 990, 350, 354, 638, 5, 27, 963, 2, 413, 605, 380, 398, 232, 547, 742, 713, 75, 616, 271, 452, 228, 979, 848, 838, 386, 251, 990, 215, 437, 444, 381, 681, 47, 393, 746, 474, 254, 267, 396, 609, 161, 963, 957, 867, 346, 634, 480, 483, 62, 315, 259, 484, 923, 327, 667, 673, 826, 74, 771, 65, 335, 308, 409, 518, 236, 795, 20, 879, 551, 928, 843, 237, 577, 995, 1009, 868, 305, 611, 405, 477, 108, 824, 876, 960, 467, 893, 205, 188, 331, 629, 150, 641, 549, 892, 635, 177, 1009, 885, 464, 1022, 796, 953, 760, 11, 781, 700, 772, 708, 315, 843, 279, 10, 17, 54, 456, 49, 904, 4, 542, 826, 871, 592, 223, 643, 587, 864, 523, 138, 934, 371, 954, 479, 337, 390, 171, 531, 377, 520, 988, 599, 452, 966, 502, 160, 215, 673, 819, 641, 340, 682, 220, 435, 386, 502, 124, 11, 966, 430, 960, 874, 835, 888, 692, 762, 369, 741, 205, 94, 49, 786, 451, 611, 757, 948, 641, 508, 654, 534, 129, 792, 968, 885, 518, 322, 630, 888, 670, 674, 259, 793, 283, 506, 447, 981, 629, 624, 407, 625, 53, 237, 238, 269, 944, 91, 161, 402, 203, 1021, 227, 554, 209, 771, 959, 880, 623, 352, 568, 148, 216, 433, 954, 130, 810, 670, 881, 616, 610, 818, 383, 955, 85, 472, 113, 385, 662, 40, 376, 361, 410, 1017, 333, 247, 934, 289, 55, 474, 367, 821, 455, 654, 686, 908, 491, 673, 357, 154, 281, 593, 567, 812, 269, 955, 57, 511, 58, 778, 113, 644, 952, 987, 836, 364, 948, 95, 866, 430, 212, 331, 402, 899, 285, 861, 34, 300, 20, 693, 558, 1021, 289, 335, 630, 833, 575, 354, 447, 85, 719, 349, 377, 928, 451, 75, 907, 399, 8, 197, 167, 583, 98, 22, 912, 429, 108, 177, 930, 743, 947, 352, 407, 15, 703, 128, 981, 681, 955, 148, 37, 388, 882, 544, 616, 250, 757, 783, 163, 748, 775, 411, 405, 830, 375, 369, 667, 237, 320, 791, 430, 446, 320, 689, 1004, 801, 59, 375, 904, 929, 793, 276, 15, 251, 935, 742, 772, 195, 870, 958, 440, 674, 739, 780, 680, 342, 693, 913, 465, 754, 757, 111, 698, 737, 627, 496, 193, 71, 817, 44, 703, 733, 521, 118, 257, 446, 654, 271, 129, 781, 196, 110, 516, 341, 376, 710, 782, 187, 820, 610, 78, 806, 845, 1004, 683, 248, 1016, 22, 693, 59, 223, 860, 605, 55, 881, 355, 902, 305, 403, 327, 859, 735, 644, 579, 955, 738, 349, 637, 39, 410, 391, 188, 258, 98, 923, 899, 680, 176, 988, 64, 655, 74, 194, 272, 125, 604, 374, 790, 415, 867, 941, 441, 969, 518, 561, 644, 784, 514, 465, 160, 662, 864, 785, 600, 906, 673, 308, 327, 454, 651, 77, 755, 406, 518, 804, 389, 322, 566, 182, 1012, 215, 894, 538, 29, 967, 861, 704, 855, 873, 814, 343, 853, 201, 106, 433, 474, 418, 476, 995, 444, 337, 599, 253, 561, 312, 739, 941, 861, 1014, 201, 549, 277, 602, 440, 176, 119, 961, 784, 737, 772, 854, 858, 335, 284, 748, 516, 947, 938, 958, 907, 472, 296, 666, 432, 69, 866, 820, 195, 722, 260, 752, 483, 80, 651, 667, 341, 770, 871, 910, 883, 477, 467, 734, 766, 948, 193, 110, 170, 578, 944, 251, 226, 494, 37, 1019, 154, 259, 598, 620, 464, 941, 554, 434, 746, 821, 412, 480, 613, 94, 763, 638, 712, 10, 118, 791, 892, 602, 894, 957, 818, 438, 54, 785, 835, 547, 683, 223, 747, 728, 175, 319, 982, 1, 757, 199, 714, 703, 308, 226, 562, 419, 789, 570, 985, 177, 495, 804, 538, 662, 193, 424, 114, 860, 1022, 371, 116, 190, 818, 298, 853, 778, 271, 243, 593, 490, 335, 785, 343, 334, 695, 781, 406, 619, 541, 862, 758, 498, 222, 633, 780, 1004, 1002, 751, 894, 60, 236, 525, 833, 999, 269, 16, 68, 798, 600, 161, 40, 150, 733, 902, 1003, 247, 77, 754, 578, 698, 670, 216, 859, 858, 309, 151, 969, 44, 708, 196, 894, 825, 170, 334, 553, 394, 575, 846, 849, 103, 339, 966, 772, 982, 972, 763, 733, 429, 357, 950, 845, 985, 317, 791, 223, 178, 687, 865, 116, 1021, 836, 296, 729, 754, 275, 155, 951, 120, 354, 441, 243, 623, 633, 326, 209, 425, 704, 605, 814, 500, 30, 871, 713, 1015, 256, 640, 29, 474, 741, 641, 193, 738, 296, 750, 74, 459, 776, 810, 339, 822, 708, 627, 882, 386, 486, 693, 873, 637, 837, 948, 652, 58, 418, 695, 850, 512, 567, 475, 781, 678, 491, 483, 1000, 423, 60, 918, 377, 148, 549, 619, 89, 592, 409, 625, 860, 755, 892, 880, 640, 203, 725, 379, 111, 390, 501, 447, 118, 635, 750, 605, 167, 595, 245, 930, 389, 149, 552, 733, 30, 684, 502, 743, 249, 431, 969, 77, 299, 232, 277, 373, 206, 745, 678, 356, 59, 446, 447, 409, 27, 634, 47, 240, 577, 217, 525, 310, 858, 550, 714, 595, 219, 517, 301, 592, 5, 319, 222, 752, 707, 682, 629, 959, 849, 220, 992, 392, 386, 429, 142, 258, 469, 542, 88, 301, 713, 507, 525, 156, 933, 883, 236, 479, 514, 374, 892, 427, 683, 571, 538, 88, 32, 37, 136, 302, 395, 618, 775, 259, 322, 257, 80, 432, 300, 651, 525, 788, 187, 997, 97, 668, 494, 340, 154, 453, 595, 331, 819, 392, 707, 575]
        sage: polynomial = truth_table_to_univariate(10, tt); polynomial
        (a^7 + a^5 + a)*x^768 + (a^8 + a^6 + a^4 + a^3 + 1)*x^640 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a)*x^576 + (a^7 + a^6 + a^4 + a^2 + a + 1)*x^544 + (a^9 + a^8 + a^6 + a^5 + a^4 + a^2 + 1)*x^528 + (a^9 + a^8 + a^7 + a^5 + a^4 + a^3 + a^2 + 1)*x^520 + (a^8 + a^6 + a^3 + a)*x^516 + (a^6 + a^5 + a^4 + a^2 + a)*x^514 + (a^7 + a^6 + a^5 + a^4 + a^3 + a^2)*x^513 + (a^9 + a^7 + a^6 + a^3 + a^2 + a)*x^512 + (a^5 + a^4 + 1)*x^384 + (a^9 + a^7 + a^5 + a^4 + a + 1)*x^320 + (a^9 + a^7 + a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^288 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a)*x^272 + (a^5 + a^4 + a)*x^264 + (a^8 + a^5 + a^3 + 1)*x^260 + (a^9 + a^8 + a^7 + a^4 + a^3 + 1)*x^258 + (a^7 + a^6 + a^5 + a^4 + a^2 + a + 1)*x^257 + (a^5 + a^4 + 1)*x^256 + (a^9 + a^8 + a^4 + 1)*x^192 + (a^9 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^160 + (a^8 + a^6 + a^4 + a^2 + 1)*x^144 + (a^9 + a^7 + a^6)*x^136 + (a^8 + a^7 + a^6 + a^3 + a + 1)*x^132 + (a^9 + a^8 + a^7 + a^5 + a^4 + a^3 + a^2)*x^130 + (a^6 + a^5 + a^4 + a^3 + a)*x^129 + (a^8 + a^5 + a^4 + a + 1)*x^128 + (a^9 + a^8 + a^6 + a^5 + 1)*x^96 + (a^8 + a^7 + a^5 + a + 1)*x^80 + (a^9 + a^8 + a^4 + a^3 + 1)*x^72 + (a^6 + a^5 + a^4 + a^3 + 1)*x^68 + (a^9 + a^8 + a^7 + a^5 + a^4 + a^2 + 1)*x^66 + (a^9 + a^8 + a^6 + a^5 + a^4)*x^65 + (a^9 + a^7 + a^6 + 1)*x^64 + (a^8 + a^7 + a^6 + a^5)*x^48 + (a^9 + a^7 + a^6 + a^3)*x^40 + (a^9 + a^6 + a^5 + a)*x^36 + (a^8 + a^6 + 1)*x^34 + (a^8 + a^7 + a^6 + a^5 + a^3 + a)*x^33 + (a^6 + a^4 + a^3 + 1)*x^32 + (a^7 + a^6 + a^4 + a^2 + a + 1)*x^24 + (a^8 + a^6 + a^5 + a^3 + 1)*x^20 + (a^7 + a^6 + a^4 + a^3)*x^18 + (a^9 + a^2 + a)*x^17 + (a^3 + a + 1)*x^16 + (a^8 + a^6 + a^4 + a^3 + a^2 + 1)*x^12 + (a^9 + a^8 + a^7 + a^6 + a^5 + a^4 + a)*x^10 + (a^8 + a^5 + a^4 + a^3 + a^2 + 1)*x^9 + (a^8 + a^4 + a^3 + a^2 + a)*x^8 + (a^7 + a^4)*x^6 + (a^7 + a^6 + a^5 + a^2 + 1)*x^5 + (a^9 + a^7 + a^6 + a^4 + a^2 + 1)*x^4 + (a^5 + a^4 + a^2)*x^3 + (a^7 + a^6 + a^5 + a^3 + a)*x^2 + (a^8 + 1)*x    
    """
    F = GF(2**n, 'a')
    return interpolate_truth_table(F, 'x', tt)


def univariate_to_matrix(n, polynomial, basis=None, output_format='univariate'):
    r"""
    Compute the Quadratic Matrix (QM) of a univariate polynomial over GF(2^n) with respect to a normal basis.
    Defined by `M[i,j] = f(b_i + b_j) + f(b_i) + f(b_j) + f(0)` for `i != j` and `M[i,i] = 0`, where `b_0, ..., b_{n-1}` is the supplied normal basis.
    
    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)
    - ``basis`` -- (optional) basis of GF(2^n) over GF(2); if None, a normal basis is used
    - ``output_format`` -- (optional) the format of the output matrix entries, either 'univariate' or 'power'; if None, defaults to 'univariate'

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import univariate_to_matrix
        sage: F.<a> = GF(2^4)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^3
        sage: M = univariate_to_matrix(4, polynomial); M
        [            0 a^3 + a^2 + a   a^2 + a + 1 a^3 + a^2 + 1]
        [a^3 + a^2 + a             0   a^3 + a + 1       a^2 + a]
        [  a^2 + a + 1   a^3 + a + 1             0       a^3 + 1]
        [a^3 + a^2 + 1       a^2 + a       a^3 + 1             0]

        sage: basis = [(a)**(i) for i in range(4)]
        sage: M = univariate_to_matrix(4, polynomial, basis, output_format='power'); M
        [   0  a^5 a^10  a^2]
        [ a^5    0  a^8 a^13]
        [a^10  a^8    0 a^11]
        [ a^2 a^13 a^11    0]

        sage: F.<a> = GF(2^6)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = a*x^24 + x^10 + x^3
        sage: M = univariate_to_matrix(6, polynomial); M
        [                        0                       a^4       a^5 + a^4 + a^2 + a                       a^2                 a^5 + a^3                   a^5 + 1]
        [                      a^4                         0                 a^5 + a^2           a^4 + a^3 + a^2       a^5 + a^4 + a^3 + 1                   a^2 + 1]
        [      a^5 + a^4 + a^2 + a                 a^5 + a^2                         0   a^4 + a^3 + a^2 + a + 1               a^4 + a + 1 a^5 + a^4 + a^3 + a^2 + 1]
        [                      a^2           a^4 + a^3 + a^2   a^4 + a^3 + a^2 + a + 1                         0         a^4 + a^2 + a + 1             a^4 + a^2 + 1]
        [                a^5 + a^3       a^5 + a^4 + a^3 + 1               a^4 + a + 1         a^4 + a^2 + a + 1                         0       a^4 + a^3 + a^2 + a]
        [                  a^5 + 1                   a^2 + 1 a^5 + a^4 + a^3 + a^2 + 1             a^4 + a^2 + 1       a^4 + a^3 + a^2 + a                         0]
    """
    F = GF(2**n, 'a')

    if basis is None:
        basis = [(F.gen()**3)**(2**i) for i in range(n)]
        
    def _poly(i, j):
        return polynomial(i + j) + polynomial(i) + polynomial(j) + polynomial(0)
    
    if output_format == 'power':
        a = F.gen()
        a_sym = SR.var('a')

        def _to_power(val):
            if val == 0: 
                return SR(0)
            elif val == 1:
                return SR(1)
            return a_sym**val.log(a)
        
        return Matrix(SR, [[_to_power(_poly(bi, bj)) for bj in basis] for bi in basis])

    return Matrix(F, [[_poly(bi, bj) for bj in basis] for bi in basis])



def matrix_to_univariate(n, M, basis=None):
    r"""
    Recover the univariate polynomial representation of a quadratic function from its Quadratic Matrix (QM) with respect to a normal basis.
    Defined by computing the coefficient matrix `C` via `C = (M_b^{-1})^T * M * M_b^{-1}`, where `M_b[i,j] = b_i^{2^j}`, and returning `f(x) = sum_{i<j} C[i,j] * x^{2^i + 2^j}`.
    
    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``M`` -- a matrix over GF(2^n)
    - ``basis`` -- (optional) basis of GF(2^n) over GF(2); if None, a normal basis is used
    
    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import matrix_to_univariate
        sage: F.<a> = GF(2^4)
        sage: M = Matrix(F, [[0, a^3 + a^2 + a, a^2 + a + 1, a^3 + a^2 + 1], [a^3 + a^2 + a, 0, a^3 + a + 1, a^2 + a], [ a^2 + a + 1, a^3 + a + 1, 0, a^3 + 1], [a^3 + a^2 + 1, a^2 + a, a^3 + 1, 0]])
        sage: polynomial = matrix_to_univariate(4, M); polynomial
        x^3

        sage: F.<a> = GF(2^6)
        sage: M = Matrix(F, [[0, a^4, a^5 + a^4 + a^2 + a, a^2, a^5 + a^3, a^5 + 1], [a^4, 0, a^5 + a^2, a^4 + a^3 + a^2, a^5 + a^4 + a^3 + 1, a^2 + 1], [a^5 + a^4 + a^2 + a, a^5 + a^2, 0, a^4 + a^3 + a^2 + a + 1, a^4 + a + 1, a^5 + a^4 + a^3 + a^2 + 1], [a^2, a^4 + a^3 + a^2, a^4 + a^3 + a^2 + a + 1, 0, a^4 + a^2 + a + 1, a^4 + a^2 + 1], [a^5 + a^3, a^5 + a^4 + a^3 + 1, a^4 + a + 1, a^4 + a^2 + a + 1, 0, a^4 + a^3 + a^2 + a], [a^5 + 1, a^2 + 1, a^5 + a^4 + a^3 + a^2 + 1, a^4 + a^2 + 1, a^4 + a^3 + a^2 + a, 0]])
        sage: polynomial = matrix_to_univariate(6, M); polynomial
        a*x^24 + x^10 + x^3
    """
    F = GF(2**n, 'a')
    R = PolynomialRing(F, 'x')
    x = R.gen()

    if M.base_ring() != F:
        raise ValueError("Matrix base ring must match the field F")
    
    if basis is None:
        basis = [(F.gen()**3)**(2**i) for i in range(n)]
        
    M_alpha = Matrix(F, [[b**(2**j) for b in basis] for j in range(n)])
    CF = M_alpha.inverse().transpose() * M * M_alpha.inverse()
    
    polynomial = R.zero()
    for i in range(n):
        for j in range(i):
            if CF[i,j] != 0:
                polynomial += CF[i,j] * x**(2**i + 2**j)

    return polynomial


def univariate_to_sequence(n, polynomial, basis=None):
    r"""
    Compute the sequence of the Quadratic Matrix (QM) of a univariate polynomial over GF(2^n) with respect to a normal basis.
    Defined by constructing the QM `M` with respect to a normal basis and reading off the upper-triangular entries row by row as integers, `[M[i,j].to_integer() for 0 <= i < j <= n-1]`.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``polynomial`` -- a univariate polynomial over GF(2^n)
    - ``basis`` -- (optional) basis of GF(2^n) over GF(2); if None, a normal basis is used

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import univariate_to_sequence
        sage: F.<a> = GF(2^4)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = x^3
        sage: sequence = univariate_to_sequence(4, polynomial); sequence
        [14, 7, 13, 11, 6, 9]

        sage: F.<a> = GF(2^6)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = a*x^24 + x^10 + x^3
        sage: sequence = univariate_to_sequence(6, polynomial); sequence
        [16, 54, 4, 40, 33, 36, 28, 57, 5, 31, 19, 61, 23, 21, 30]

        sage: F.<a> = GF(2^8)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = (a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^192 + (a^6 + a^5 + a^2 + 1)*x^160 + (a^7 + a^6 + a^4 + 1)*x^144 + (a^6 + a^4 + a^3 + 1)*x^136 + (a^7 + a^5 + a + 1)*x^132 + (a^7 + a^6 + a^4 + a^3 + a^2 + a)*x^130 + (a^6 + a^4 + a^2 + a + 1)*x^129 + (a^7 + a^6 + a^4 + a^2)*x^96 + (a^5 + a^3 + a + 1)*x^80 + (a^6 + a^2)*x^72 + (a^7 + a^5 + a)*x^68 + a^7*x^66 + (a^7 + a^6 + a^5 + a^2 + a)*x^65 + (a^7 + a^4 + a^3 + a^2 + a + 1)*x^48 + (a^7 + a^3 + a^2 + a)*x^40 + (a^3 + 1)*x^36 + (a^7 + a^6 + a^5)*x^34 + (a^6 + a)*x^33 + (a^6 + a^4 + a)*x^24 + (a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^20 + (a^7 + a^6 + a^5 + a^3 + a + 1)*x^18 + (a^7 + a^4 + a + 1)*x^17 + (a^7 + a^5 + a^2 + 1)*x^12 + (a^6 + a^5 + a^4 + a^3 + a + 1)*x^10 + (a^4 + a^2 + a)*x^9 + a*x^6 + (a^6 + a^5 + a^4 + a^3 + a^2 + a)*x^5 + (a^7 + a^4 + a^3)*x^3
        sage: basis = [(a)^(i) for i in range(8)]
        sage: sequence = univariate_to_sequence(8, polynomial, basis); sequence
        [1, 6, 204, 20, 142, 72, 85, 204, 202, 154, 20, 85, 29, 8, 48, 30, 160, 61, 30, 46, 157, 160, 34, 39, 209, 209, 246, 175]
    """
    F = GF(2**n, 'a')

    if basis is None:
        basis = [(F.gen()**3)**(2**i) for i in range(n)]

    M = univariate_to_matrix(n, polynomial, basis)
    return [M[i, j].to_integer() for i in range(n) for j in range(i + 1, n)]


def sequence_to_univariate(n, sequence, basis=None):
    r"""
    Recover the univariate polynomial representation of a sequence of the Quadratic Matrix (QM) with respect to a normal basis.
    Defined by filling the symmetric matrix `M` with `M[i,j] = M[j,i] = s_k` for the k-th sequence element (in row-major upper-triangular order), then applying `matrix_to_polynomial` to recover `f`.
    
    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n)
    - ``sequence`` -- a list of integers representing the upper-triangular entries of the QM, of length `n*(n-1)/2`
    - ``basis`` -- (optional) basis of GF(2^n) over GF(2); if None, a normal basis is used

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import sequence_to_univariate
        sage: sequence = [14, 7, 13, 11, 6, 9]
        sage: polynomial = sequence_to_univariate(4, sequence); polynomial
        x^3

        sage: sequence = [16, 54, 4, 40, 33, 36, 28, 57, 5, 31, 19, 61, 23, 21, 30]
        sage: polynomial = sequence_to_univariate(6, sequence); polynomial
        a*x^24 + x^10 + x^3

        sage: F.<a> = GF(2^8)
        sage: sequence = [1, 6, 204, 20, 142, 72, 85, 204, 202, 154, 20, 85, 29, 8, 48, 30, 160, 61, 30, 46, 157, 160, 34, 39, 209, 209, 246, 175]
        sage: basis = [(a)^(i) for i in range(8)]
        sage: polynomial = sequence_to_univariate(8, sequence, basis); polynomial
        (a^6 + a^5 + a^4 + a^3 + a^2 + a + 1)*x^192 + (a^6 + a^5 + a^2 + 1)*x^160 + (a^7 + a^6 + a^4 + 1)*x^144 + (a^6 + a^4 + a^3 + 1)*x^136 + (a^7 + a^5 + a + 1)*x^132 + (a^7 + a^6 + a^4 + a^3 + a^2 + a)*x^130 + (a^6 + a^4 + a^2 + a + 1)*x^129 + (a^7 + a^6 + a^4 + a^2)*x^96 + (a^5 + a^3 + a + 1)*x^80 + (a^6 + a^2)*x^72 + (a^7 + a^5 + a)*x^68 + a^7*x^66 + (a^7 + a^6 + a^5 + a^2 + a)*x^65 + (a^7 + a^4 + a^3 + a^2 + a + 1)*x^48 + (a^7 + a^3 + a^2 + a)*x^40 + (a^3 + 1)*x^36 + (a^7 + a^6 + a^5)*x^34 + (a^6 + a)*x^33 + (a^6 + a^4 + a)*x^24 + (a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^20 + (a^7 + a^6 + a^5 + a^3 + a + 1)*x^18 + (a^7 + a^4 + a + 1)*x^17 + (a^7 + a^5 + a^2 + 1)*x^12 + (a^6 + a^5 + a^4 + a^3 + a + 1)*x^10 + (a^4 + a^2 + a)*x^9 + a*x^6 + (a^6 + a^5 + a^4 + a^3 + a^2 + a)*x^5 + (a^7 + a^4 + a^3)*x^3
    """
    F = GF(2**n, 'a')
    
    if basis is None:
        basis = [(F.gen()**3)**(2**i) for i in range(n)]
    
    seq_iter = (F.from_integer(s) for s in sequence)
    M = Matrix(F, n)
    for i in range(n):
        for j in range(i + 1, n):
            M[i,j] =  M[j,i] = next(seq_iter)

    return matrix_to_univariate(n, M, basis)


def univariate_to_bivariate(n, polynomial, basis=None):
    r"""
    Convert a polynomial from its univariate representation over GF(2^(2m)) to its bivariate representation over GF(2^m).
    Defined by `F(x, y) = (f(x, y), g(x, y))`.

    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^(2m)); must be even
    - ``polynomial`` -- a univariate polynomial over GF(2^(2m))
    - ``basis`` -- (optional) basis (b0, b1) of GF(2^(2m)) over GF(2^m); if None, a canonical basis (1, a) is used, where a is the generator of GF(2^(2m))

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import univariate_to_bivariate
        sage: F.<a> = GF(2^6)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = (a^4 + a^3 + 1)*x^17 + (a + 1)*x^8 + (a^4 + a^3)*x^3 + x
        sage: f, g = univariate_to_bivariate(6, polynomial); f, g
        (x^3 + x*y^2, x^2*y + y^3 + x)

        sage: polynomial = (a^4 + a^3 + 1)*x^16 + (a^4 + a^2 + a)*x^9 + (a^3 + a^2 + a)*x^2 + a
        sage: f, g = univariate_to_bivariate(6, polynomial); f, g
        (x^2 + x*y + y^2, y^2 + 1)
    """
    if n % 2 != 0:
        raise ValueError("n must be even to convert to bivariate representation")
    
    m = n // 2
    q = 2**m

    F = GF(2**n, 'a')
    a = F.gen()
    F_sub, hom = F.subfield(m, 'b', map=True)
    hom_inv = hom.section()

    if basis is None:
        basis = [F(1), a]
    b0, b1 = basis[0], basis[1]

    univariate_tt = univariate_to_truth_table(n, polynomial)
    delta = b0**q * b1 + b0 * b1**q

    # Build the two bivariate truth tables, indexed by integer codes of (x, y)
    f_tt = [[0] * q for _ in range(q)]
    g_tt = [[0] * q for _ in range(q)]
    for i in range(q):
        FX = hom(F_sub.from_integer(i))
        for j in range(q):
            FY = hom(F_sub.from_integer(j))
            X = b0 * FX + b1 * FY
            c = F.from_integer(univariate_tt[X.to_integer()])
            # Solve the system of equations
            c1 = (b0**q * c + b0 * c**q) / delta
            c0 = (c + c1 * b1) / b0
            f_tt[i][j] = hom_inv(c0).to_integer()
            g_tt[i][j] = hom_inv(c1).to_integer()

    R = PolynomialRing(F_sub, ['x', 'y'])
    x, y = R.gens()

    # Construct the bivariate polynomials via Lagrange interpolation on the bivariate truth tables, using the isomorphism hom to read off coefficients
    def interpolate(tt):
        layer_y = [interpolate_truth_table(F_sub, 'y', tt[i]) for i in range(q)]

        result = R(0)
        for dy in range(q):
            coeff_tt = [layer_y[i][dy].to_integer() for i in range(q)]
            poly_in_x = interpolate_truth_table(F_sub, 'x', coeff_tt)
            for dx, c in enumerate(poly_in_x.list()):
                if c != 0:
                    result += c * x**dx * y**dy
        return result

    return interpolate(f_tt), interpolate(g_tt)


def bivariate_to_univariate(m, f, g, basis=None):
    r"""
    Convert a pair of polynomials from their bivariate representation over GF(2^m) to their univariate representation over GF(2^(2m)).
    Defined by `F(X) = b0 * f(x) + b1*g(x)`, where `b0, b1` are elements of the basis of GF(2^(2m)) over GF(2^m) and `f(x), g(x)` are the bivariate polynomials evaluated at the coordinates of `x` under the isomorphism between GF(2^(2m)) and GF(2^m)^2.

    INPUT:

    - ``m`` -- the degree of the finite field extension GF(2^m)
    - ``f`` -- a bivariate polynomial over GF(2^m)
    - ``g`` -- a bivariate polynomial over GF(2^m)
    - ``basis`` -- (optional) basis (b0, b1) of GF(2^(2m)) over GF(2^m); if None, a canonical basis (1, a) is used, where a is the generator of GF(2^(2m))

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import bivariate_to_univariate
        sage: F.<b> = GF(2^3)
        sage: R.<x, y> = PolynomialRing(F)
        sage: f = x^3 + x*y^2
        sage: g = x^2*y + y^3 + x
        sage: polynomial = bivariate_to_univariate(3, f, g); polynomial
        (a^4 + a^3 + 1)*x^17 + (a + 1)*x^8 + (a^4 + a^3)*x^3 + x

        sage: f = x^2 + x*y + y^2
        sage: g = y^2 + 1
        sage: polynomial = bivariate_to_univariate(3, f, g); polynomial
        (a^4 + a^3 + 1)*x^16 + (a^4 + a^2 + a)*x^9 + (a^3 + a^2 + a)*x^2 + a
    """
    n = 2 * m
    q = 2**m
    
    F = GF(2**n, 'a')
    a = F.gen()
    F_sub, hom = F.subfield(m, 'b', map=True)

    if basis is None:
        basis = [F(1), a]
    b0, b1 = basis[0], basis[1]
    
    # Iterate over all pairs (x, y) in the small field, evaluate f, g, and fill in the univariate truth table via the isomorphism hom
    univariate_tt = [0] * (q ** 2)
    for i in range(q):
        x_sub = F_sub.from_integer(i)
        for j in range(q):
            y_sub = F_sub.from_integer(j)
            X = b0 * hom(x_sub) + b1 * hom(y_sub)
            f_val = hom(f(x_sub, y_sub))
            g_val = hom(g(x_sub, y_sub))
            FX = b0 * f_val + b1 * g_val
            univariate_tt[X.to_integer()] = FX.to_integer()
    
    return interpolate_truth_table(F, 'x', univariate_tt)


def univariate_to_trivariate(n, polynomial, basis=None):
    r"""
    Convert a polynomial from its univariate representation over GF(2^(3m)) to its trivariate representation over GF(2^m).
    Defined by `F(x, y, z) = (f(x, y, z), g(x, y, z), h(x, y, z))`.
    
    INPUT:

    - ``n`` -- the degree of the finite field extension GF(2^n); must be divisible by 3
    - ``polynomial`` -- a univariate polynomial over GF(2^n)
    - ``basis`` -- (optional) basis (b0, b1, b2) of GF(2^(3m)) over GF(2^m); if None, a canonical basis (1, a, a^2) is used, where a is the generator of GF(2^(3m))

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import univariate_to_trivariate
        sage: F.<a> = GF(2^9)
        sage: R.<x> = PolynomialRing(F)
        sage: polynomial = (a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + a)*x^192 + (a^8 + a^7 + a^6 + a^3 + a)*x^136 + (a^7 + a^6 + a^3 + a^2 + a)*x^129 + (a^6 + a^4 + a^3 + 1)*x^80 + (a^5 + a + 1)*x^66 + (a^8 + a^6 + a^4 + a^2 + a)*x^24 + (a^7 + a^5 + a^3 + a^2)*x^17 + (a^5 + a^3 + a^2 + a + 1)*x^10 + (a^8 + a^6 + a^2)*x^3
        sage: f, g, h = univariate_to_trivariate(9, polynomial); f, g, h
        (x^3 + x^2*z + y*z^2, y^3 + x^2*z, x*y^2 + y^2*z + z^3)

        sage: polynomial = (a^7 + a^5 + a^3 + 1)*x^192 + (a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^136 + (a^7 + a^3 + 1)*x^129 + (a^6 + a^4 + 1)*x^80 + (a^8 + a^5 + a^3 + a)*x^66 + (a^8 + a^6 + a^5 + a + 1)*x^24 + (a^8 + a^6 + a^5 + a^2 + a + 1)*x^17 + (a^6 + a^5 + 1)*x^10 + (a^8 + a^7 + a^6 + a^5 + a)*x^3
        sage: f, g, h = univariate_to_trivariate(9, polynomial); f, g, h
        (x^3 + x*y^2 + y*z^2, x*y^2 + z^3, y^3 + x^2*z + y^2*z)
    """
    if n % 3 != 0:
        raise ValueError("n must be divisible by 3 to convert to trivariate representation")
    
    m = n // 3
    q = 2**m

    F = GF(2**n, 'a')
    a = F.gen()
    F_sub, hom = F.subfield(m, 'b', map=True)
    hom_inv = hom.section()

    if basis is None:
        basis = [F(1), a, a**2]
    b0, b1, b2 = basis[0], basis[1], basis[2]

    univariate_tt = univariate_to_truth_table(n, polynomial)
    M_inv = Matrix(F, [[b0,         b1,         b2        ],
                       [b0**q,      b1**q,      b2**q     ],
                       [b0**(q**2), b1**(q**2), b2**(q**2)]]).inverse()

    # Build the three trivariate truth tables, indexed by integer codes of (x, y, z)
    f_tt = [[[0]*q for _ in range(q)] for _ in range(q)]
    g_tt = [[[0]*q for _ in range(q)] for _ in range(q)]
    h_tt = [[[0]*q for _ in range(q)] for _ in range(q)]
    for i in range(q):
        FX = hom(F_sub.from_integer(i))
        for j in range(q):
            FY = hom(F_sub.from_integer(j))
            for k in range(q):
                FZ = hom(F_sub.from_integer(k))
                X = b0 * FX + b1 * FY + b2 * FZ
                c = F.from_integer(univariate_tt[X.to_integer()])
                # Solve the system of equations
                c0, c1, c2 = M_inv * vector(F, [c, c**q, c**(q**2)])
                f_tt[i][j][k] = hom_inv(c0).to_integer()
                g_tt[i][j][k] = hom_inv(c1).to_integer()
                h_tt[i][j][k] = hom_inv(c2).to_integer()

    R = PolynomialRing(F_sub, ['x', 'y', 'z'])
    x, y, z = R.gens()

    # Construct the trivariate polynomials via Lagrange interpolation on the trivariate truth tables, using the isomorphism hom to read off coefficients
    def interpolate(tt):
        layer_z = [[interpolate_truth_table(F_sub, 'z', tt[i][j]) for j in range(q)] for i in range(q)]
        layer_yz = [[interpolate_truth_table(F_sub, 'y', [layer_z[i][j][dz].to_integer() for j in range(q)]) for dz in range(q)] for i in range(q)]
        
        result = R(0)
        for dz in range(q):
            for dy in range(q):
                coeff_tt = [layer_yz[i][dz][dy].to_integer() for i in range(q)]
                poly_in_x = interpolate_truth_table(F_sub, 'x', coeff_tt)
                for dx, c in enumerate(poly_in_x.list()):
                    if c != 0:
                        result += c * x**dx * y**dy * z**dz
        return result

    return interpolate(f_tt), interpolate(g_tt), interpolate(h_tt)

    
def trivariate_to_univariate(m, f, g, h, basis=None):
    r"""
    Convert a triplet of polynomials from their trivariate representation over GF(2^m) to their univariate representation over GF(2^n).
    Defined by `F(X) = b0*f(x) + b1*g(x) + b2*h(x)`, where `b0, b1, b2` are elements of the basis of GF(2^n) over GF(2^m) and `f(x), g(x), h(x)` are the trivariate polynomials evaluated at the coordinates of `x` under the isomorphism between GF(2^n) and GF(2^m)^3.

    INPUT:

    - ``m`` -- the degree of the finite field extension GF(2^m)
    - ``f`` -- a trivariate polynomial over GF(2^m)
    - ``g`` -- a trivariate polynomial over GF(2^m)
    - ``h`` -- a trivariate polynomial over GF(2^m)
    - ``basis`` -- (optional) basis (b0, b1, b2) of GF(2^(3m)) over GF(2^m); if None, a canonical basis (1, a, a^2) is used, where a is the generator of GF(2^(3m))

    EXAMPLES::

        sage: from sagemath_cryptographic_functions_library import trivariate_to_univariate
        sage: F.<b> = GF(2^3)
        sage: R.<x, y, z> = PolynomialRing(F)
        sage: f = x^3 + x^2*z + y*z^2
        sage: g = y^3 + x^2*z
        sage: h = x*y^2 + y^2*z + z^3
        sage: polynomial = trivariate_to_univariate(3, f, g, h); polynomial
        (a^8 + a^7 + a^6 + a^5 + a^3 + a^2 + a)*x^192 + (a^8 + a^7 + a^6 + a^3 + a)*x^136 + (a^7 + a^6 + a^3 + a^2 + a)*x^129 + (a^6 + a^4 + a^3 + 1)*x^80 + (a^5 + a + 1)*x^66 + (a^8 + a^6 + a^4 + a^2 + a)*x^24 + (a^7 + a^5 + a^3 + a^2)*x^17 + (a^5 + a^3 + a^2 + a + 1)*x^10 + (a^8 + a^6 + a^2)*x^3
        
        sage: f = x^3 + x*y^2 + y*z^2
        sage: g = x*y^2 + z^3
        sage: h = y^3 + x^2*z + y^2*z
        sage: polynomial = trivariate_to_univariate(3, f, g, h); polynomial
        (a^7 + a^5 + a^3 + 1)*x^192 + (a^7 + a^6 + a^4 + a^3 + a^2 + 1)*x^136 + (a^7 + a^3 + 1)*x^129 + (a^6 + a^4 + 1)*x^80 + (a^8 + a^5 + a^3 + a)*x^66 + (a^8 + a^6 + a^5 + a + 1)*x^24 + (a^8 + a^6 + a^5 + a^2 + a + 1)*x^17 + (a^6 + a^5 + 1)*x^10 + (a^8 + a^7 + a^6 + a^5 + a)*x^3
    """
    n = 3 * m
    q = 2**m
    
    F = GF(2**n, 'a')
    a = F.gen()
    F_sub, hom = F.subfield(m, 'b', map=True)

    if basis is None:
        basis = [F(1), a, a**2]
    b0, b1, b2 = basis[0], basis[1], basis[2]
    
    # Iterate over all triples (x, y, z) in the small field, evaluate f, g, h, and fill in the univariate truth table via the isomorphism hom.
    univariate_tt = [0] * (q ** 3)
    for i in range(q):
        x_sub = F_sub.from_integer(i)
        for j in range(q):
            y_sub = F_sub.from_integer(j)
            for k in range(q):
                z_sub = F_sub.from_integer(k)
                X = b0 * hom(x_sub) + b1 * hom(y_sub) + b2 * hom(z_sub)
                f_val = hom(f(x_sub, y_sub, z_sub))
                g_val = hom(g(x_sub, y_sub, z_sub))
                h_val = hom(h(x_sub, y_sub, z_sub))
                FX = b0 * f_val + b1 * g_val + b2 * h_val
                univariate_tt[X.to_integer()] = FX.to_integer()
    
    return interpolate_truth_table(F, 'x', univariate_tt)
