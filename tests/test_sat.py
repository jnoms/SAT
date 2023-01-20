from sat import __version__


def test_version():
    assert __version__ == "0.1.0"


def test_dicts_are_ordered():
    """
    This is a major assumption of some functions. So, need to catch cases where folks
    are trying to use a python version that doesn't use ordered dictionaries.
    """

    test_dict = dict()
    test_dict[0] = "one"
    test_dict[1] = "two"
    test_dict[2] = "three"

    key_list = list(test_dict.keys())

    assert key_list[0] == 0
    assert key_list[1] == 1
    assert key_list[2] == 2
