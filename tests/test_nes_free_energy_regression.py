import pathlib
from contextlib import redirect_stdout
from io import StringIO

from water_nes.analysis.nes_free_energy import calculate_nes_free_energy


def test_pmx_calc_sep_regression(data_regression, file_regression):
    # The test input is located in a folder called nes_input_files/ in the same folder
    # as the test file is located. We're then using the relative path compared to the
    # current work directory, such that the path is going to look the same on every
    # system, as long as the test is run from the same directory (typically, the
    # root directory of the repo).
    test_directory = pathlib.PurePath(__file__).parent
    input_directory = test_directory.joinpath(
        test_directory, "nes_input_files"
    ).relative_to(pathlib.Path.cwd())

    # Create lists of input files. The transition between stages A and B happens in
    # two steps, by first turning off the Coulomb interactions, then the vdW
    # interactions, or vice-versa.
    xvg_forward_coulomb = [
        input_directory.joinpath(f"transition_A2B_coul_{n}.xvg") for n in range(1, 11)
    ]
    xvg_forward_vdw = [
        input_directory.joinpath(f"transition_A2B_vdw_{n}.xvg") for n in range(1, 11)
    ]
    xvg_backward_coulomb = [
        input_directory.joinpath(f"transition_B2A_coul_{n}.xvg") for n in range(1, 11)
    ]
    xvg_backward_vdw = [
        input_directory.joinpath(f"transition_B2A_vdw_{n}.xvg") for n in range(1, 11)
    ]

    # Redirect stdout into string which we can test
    test_output = StringIO()
    with redirect_stdout(test_output):
        free_energy_estimate = calculate_nes_free_energy(
            xvg_files_forward_transition=[xvg_forward_coulomb, xvg_forward_vdw],
            xvg_files_backward_transition=[xvg_backward_vdw, xvg_backward_coulomb],
            temperature=298.15,
            output_units="kcal/mol",
            bootstrapping_repeats=0,
        )
        print(round(free_energy_estimate, 4))

    # Test return value, using a reasonable precision
    # (data_regression fixture requires dict input)
    data_regression.check(round(free_energy_estimate, 4)._asdict())
    # Test printed output
    file_regression.check(contents=test_output.getvalue())
