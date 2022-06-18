import pathlib

import pytest

from water_nes.analysis.nes_free_energy import (
    calculate_nes_free_energy,
    convert_energy,
    find_unique_list_size,
    make_list_of_lists,
)


class TestMakeListOfLists:
    r"""
    Tests the functionality of make_list_of_lists()
    """

    @staticmethod
    def test_valid_entries_simple() -> None:
        r"""
        Test entries that are expected to be valid single lists
        """
        entries = [
            # Simple empty list
            [[], []],
            # Simple list of str
            [["a", "b", "c"], ["d", "e", "f"]],
            # Simple list of paths
            [[pathlib.Path("abc")], [pathlib.Path("def")]],
            [[pathlib.PurePath("zx")], [pathlib.PurePath("we")]],
            # Simple mixed lists
            [["a", pathlib.Path("1")], [pathlib.PurePath("ns"), 3]],
        ]

        for entry in entries:
            assert make_list_of_lists(entry[0], entry[1]) == ([entry[0]], [entry[1]])

    @staticmethod
    def test_valid_entries_double() -> None:
        r"""
        Test entries that are expected to be valid lists of lists
        """
        entries = [
            # Double empty list
            [[[]], [[]]],
            # Double list of str
            [[["a", "b"], ["c", "d"]], [["e", "f"], ["g", "h"]]],
            # Double list of paths
            [
                [[pathlib.Path("abc")], [pathlib.Path("def")]],
                [[pathlib.Path("ghi")], [pathlib.Path("jkl")]],
            ],
            [[[pathlib.PurePath("zx")]], [[pathlib.PurePath("we")]]],
            # Double mixed lists
            [[["a", pathlib.Path("1")]], [[pathlib.PurePath("ns"), 3]]],
        ]

        for entry in entries:
            assert make_list_of_lists(entry[0], entry[1]) == (entry[0], entry[1])

    @staticmethod
    def test_invalid_entries() -> None:
        r"""
        Test entries that are expected to be invalid
        """
        entries = [
            # Single / double empty list
            [[], [[]]],
            # Single / double list of str
            [[["a", "b"], ["c", "d"]], ["e", "f", "g", "h"]],
            # Other types
            [1, []],
            ["s", "b"],
            [None, None],
        ]

        for entry in entries:
            with pytest.raises(TypeError):
                make_list_of_lists(entry[0], entry[1])


class TestFindUniqueListSize:
    r"""
    Tests the functionality of find_unique_list_size()
    """

    @staticmethod
    def test_valid_input() -> None:
        r"""
        Test input that is expected to yield a unique size
        """
        entries = [
            ({"a": [[], []], "b": [[], []]}, 0),
            ({"a": [[1]], "b": [[3]]}, 1),
            ({"a": [[1, 2]], "b": [[3, 4]]}, 2),
            (
                {
                    "a": [["abc", "de", "f"], [1, 2, 3]],
                    "b": [[1, 1, 1], ["xx", "xx", "xx"]],
                },
                3,
            ),
        ]
        for entry in entries:
            assert find_unique_list_size(entry[0]) == entry[1]

    @staticmethod
    def test_invalid_input() -> None:
        r"""
        Test input that is expected to raise an error
        """
        entries = [
            {"a": [[], ["x"]], "b": [[], []]},
            {"a": [[1]], "b": [[3, 2]]},
            {"a": [[1, 2]], "b": [[3, 4, 5]]},
            {"a": [["abc", "de", "f"], [1, 2]], "b": [[1, 1, 1], ["xx", "xx", "xx"]]},
        ]
        for entry in entries:
            with pytest.raises(RuntimeError):
                find_unique_list_size(entry)


class TestConvertEnergy:
    r"""
    Tests the functionality of convert_energy()
    """

    @staticmethod
    def test_kj_mol() -> None:
        r"""
        Test that conversion to kJ / mol returns identical value
        """
        for number in [-1.634e13, -162, -1.93e-97, 0, 1.23e-102, 155.781, 1.433e23]:
            assert convert_energy(number, "kJ/mol") == number

    @staticmethod
    def test_kcal_mol() -> None:
        r"""
        Test that conversion to kcal / mol exhibits expected behavior
        """
        for number in [-1.634e13, -162, -1.93e-97, 0, 1.23e-102, 155.781, 1.433e23]:
            assert convert_energy(number, "kcal/mol") == number * 0.2390057361

    @staticmethod
    def test_unimplemented() -> None:
        r"""
        Test that unimplemented conversions throw
        """
        for units in ["kj/mol", "kJ / mol", "kcal / mol", "kT", "kt", "abc"]:
            with pytest.raises(NotImplementedError):
                convert_energy(1.0, units)


class TestCalculateNesFreeEnergy:
    r"""
    Tests some part of the functionality of calculate_nes_free_energy().
    Note that there are regression tests covering more of the functionality in
    test_nes_free_energy_regression.py.
    """

    @staticmethod
    def test_type_error() -> None:
        r"""
        Tests that type error from incompatible combination of xvg input is raised
        """
        with pytest.raises(TypeError):
            calculate_nes_free_energy(
                xvg_files_forward_transition=[["1.xvg"], ["2.xvg"]],
                xvg_files_backward_transition=["1.xvg", "2.xvg"],
                temperature=300,
                output_units="kJ/mol",
                bootstrapping_repeats=0,
            )
        with pytest.raises(TypeError):
            calculate_nes_free_energy(
                xvg_files_forward_transition=[["1.xvg", "2.xvg"]],
                xvg_files_backward_transition=[["1.xvg", "2.xvg", "3.xvg"]],
                temperature=300,
                output_units="kJ/mol",
                bootstrapping_repeats=0,
            )
