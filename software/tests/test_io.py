import pytest
import os
import yaml
from diaux.io import standardize_strains, scrape_frontmatter

def test_standardize_strains():
    """
    Test the standardize_strains function from the io module.

    The standardize_strains function takes a list of GC Database strain identifiers 
    and returns a standardized list of lists with the standardized shorthand notations,
    the strain genotype, the annotation of the lab stock, and the strain class.

    This test checks that the function returns the correct output structure for a 
    given input. It does not check the actual values returned by the function, as 
    these would depend on the specific strains input and the current state of the 
    __STRAINS__ global variable.
    """
    # Define a list of strains for testing
    strains = ['GC001', 'GC029']

    # Call the function with the test strains
    result = standardize_strains(strains)

    # Check that the result is a list of three lists
    assert isinstance(result, list)
    assert len(result) == 3
    assert all(isinstance(sublist, list) for sublist in result)

    # Check that each sublist has the same length as the input list of strains
    assert all(len(sublist) == len(strains) for sublist in result)

    # Ensure that the test fails if the strain is not in the database
    try:
        standardize_strains(['GC002'])
        assert False
    except KeyError:
        assert True


    
def test_scrape_frontmatter(tmpdir):
    """
    Test the scrape_frontmatter function from the io module.

    The scrape_frontmatter function reads the status of a given experimental dataset 
    from a README.md file in a specified directory. The status is embedded in the 
    README.md file as a YAML metadata block.

    This test checks that the function correctly reads the status from a test README.md 
    file and returns it as a dictionary. It also checks that the function correctly 
    handles cases where the status is missing or not an acceptable value.
    """
    # Create a test README.md file with a YAML metadata block
    readme_content = "---\nstatus: accepted\n---\n"
    readme_file = tmpdir.join("README.md")
    readme_file.write(readme_content)

    # Call the function with the test directory
    result = scrape_frontmatter(str(tmpdir))

    # Check that the result is a dictionary with the correct status
    assert isinstance(result, dict)
    assert result == {"status": "accepted"}

    # Check that the result is rejected.
    readme_content = "---\nstatus: rejected\n---\n"
    readme_file = tmpdir.join("README.md")
    readme_file.write(readme_content)
    result = scrape_frontmatter(str(tmpdir))
    assert isinstance(result, dict)
    assert result == {"status": "rejected"}

    # Change the status to an unacceptable value and test again
    readme_content = "---\nstatus: blah\n---\n"
    readme_file = tmpdir.join("README.md")
    readme_file.write(readme_content)
    result = scrape_frontmatter(str(tmpdir))
    assert isinstance(result, dict) 
    assert result == {}

    # Remove the status and test again
    readme_content = "---\n---\n"
    readme_file = tmpdir.join("README.md")
    readme_file.write(readme_content)
    result = scrape_frontmatter(str(tmpdir))
    assert result == {} 