import pathlib

from water_nes.analysis.endpoint_free_energy import calculate_endpoint_free_energy


def test_endpoint_free_energy_regression(data_regression):
    # The test input is located in a folder called endpoint_input_files/ in the same
    # folder the test file is located. We're then using the relative path compared to
    # the current work directory, such that the path is going to look the same on every
    # system, as long as the test is run from the same directory (typically, the
    # root directory of the repo).
    test_directory = pathlib.PurePath(__file__).parent
    input_directory = test_directory.joinpath(
        test_directory, "endpoint_input_files"
    ).relative_to(pathlib.Path.cwd())

    free_energy_estimate = calculate_endpoint_free_energy(
        file_lambda_0=input_directory.joinpath("lambda0.xvg"),
        file_lambda_1=input_directory.joinpath("lambda1.xvg"),
        start_time=1000,
        end_time=15000,
        output_units="kcal/mol",
    )

    # Test return value, using a reasonable precision
    # (data_regression fixture requires dict input)
    data_regression.check(round(free_energy_estimate, 4).as_dict())
