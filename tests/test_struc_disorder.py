from sat.scripts.struc_disorder import find_disorder, find_order


def test_count_disorder():
    plddts = [
        60,
        40,
        35,
        50,
        45,
        55,
        30,
        20,
        40,
        50,
        10,
        45,
        50,
        80,
        90,
        30,
        20,
        10,
        5,
        1,
    ]
    cutoff = 50
    n_sequential = 4

    stretches = find_disorder(plddts, cutoff, n_sequential)
    residues = [item for sublist in stretches for item in sublist]
    result = len(residues)
    assert result == 16, f"Expected 16, but got {result}"


def test_count_disorder_punctuated():
    plddts = [
        60,
        40,
        35,
        50,
        45,
        55,
        80,
        20,
        40,
        50,
        10,
        45,
        50,
        80,
        90,
        30,
        20,
        10,
        5,
        1,
    ]
    cutoff = 50
    n_sequential = 6

    stretches = find_disorder(plddts, cutoff, n_sequential)
    residues = [item for sublist in stretches for item in sublist]
    result = len(residues)
    assert result == 6, f"Expected 6, but got {result}"


def test_count_order():
    plddts = [60, 40, 35, 50, 85, 90, 95, 100, 110, 80, 85, 90, 70, 75, 80, 70, 75]
    cutoff = 80
    n_sequential = 4

    stretches = find_order(plddts, cutoff, n_sequential)
    residues = [item for sublist in stretches for item in sublist]
    result = len(residues)
    assert result == 8, f"Expected 8, but got {result}"


def test_count_order_punctuated():
    plddts = [60, 40, 35, 85, 85, 90, 95, 100, 110, 40, 85, 90, 70, 75, 80, 70, 75]
    cutoff = 80
    n_sequential = 6

    stretches = find_order(plddts, cutoff, n_sequential)
    residues = [item for sublist in stretches for item in sublist]
    result = len(residues)
    assert result == 6, f"Expected 6, but got {result}"
