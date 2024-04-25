import pandas as pd
import numpy as np 
from diaux.quant import classify_diauxic_phases, compute_pearson_correlation

def test_compute_pearson_correlation():
    """
    Test the compute_pearson_correlation function to ensure it is working correctly.
    This test function checks the output type, output shape, and a simple test case with known output.
    """
    # Test data
    data = pd.DataFrame({
        'time_hr': np.arange(100),
        'od650nm': np.random.rand(100)
    })

    # Call the function
    pearson_window=10
    result = compute_pearson_correlation(data,
                                         pearson_window=pearson_window)

    # Check the output type
    assert isinstance(result, pd.DataFrame)

    # Check the output shape
    assert result.shape[0] == data.shape[0] - pearson_window

    # Test with known output
    data = pd.DataFrame({
        'time_hr': np.arange(100),
        'od650nm': np.exp(np.linspace(1, 10, 100))
    })
    result = compute_pearson_correlation(data, savgol_filter=False)
    print(result['pearson_correlation'].max(), result['pearson_correlation'].min())
    assert np.all(np.isclose(result['pearson_correlation'].values, 1, atol=1E-4))

def test_classify_diauxic_phases():
    """
    Test the classify_diauxic_phases function to ensure it is working correctly.
    This test function checks the output type, output shape, and a simple test case with known output.
    """
    # Test data with four phases
    data = pd.DataFrame({
        'time_hr': np.arange(200),
        'pearson_correlation': np.concatenate([np.ones(50), np.zeros(50), np.full(50, 0.5), np.ones(50)]),
        'od650nm': np.random.rand(200)
    })

    # Call the function
    result = classify_diauxic_phases(data)

    # Check the output type
    assert isinstance(result, pd.DataFrame)

    # Check the output shape
    assert result.shape[0] == data.shape[0]

    # Test with known output
    assert 'phase' in result.columns
    assert set(result['phase'].unique()) == {'preshift_exponential', 'postshift_exponential', 'stall', 'transition'} 


