from sat.scripts.struc_detect_interaction import domain_chian_counts


def test_domain_chian_counts_ordered():
    assert domain_chian_counts(frozenset([0, 1, 2, 3, 4, 5]), 3) == (3, 3)


def test_domain_chian_counts_high_values():
    assert domain_chian_counts(frozenset([9, 19, 29, 39, 49, 59]), 29) == (2, 4)


def test_domain_chian_counts_unordered():
    assert domain_chian_counts(frozenset([4, 0, 5, 1, 2, 3]), 3) == (3, 3)


def test_domain_chian_counts_no_greater():
    assert domain_chian_counts(frozenset([0, 1, 2]), 3) == (3, 0)


def test_domain_chian_counts_no_smaller():
    assert domain_chian_counts(frozenset([3, 4, 5]), 3) == (0, 3)
