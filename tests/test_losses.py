"""
Tests for energy loss functions in plasma_sim.losses.

Each test verifies a single physical property of the loss functions.
"""

import numpy as np
from plasma_sim.losses import synchrotron_loss, inverse_compton_loss, bremsstrahlung_loss, coulomb_loss

"SYNCHROTRON LOSS TESTS"


def test_synchrotron_loss_is_negative():
    """Test that synchrotron loss always removes energy from the el..

    Given: an el. with a typical Lorentz factor and a typical magnetic field
    When: the synchrotron loss rate is computed
    Then: the result must be negative since el. only lose energy
    """
    assert synchrotron_loss(gamma=1e4, B=1e-9) < 0


def test_synchrotron_loss_scales_with_gamma_squared():
    """Test that synchrotron loss scales as gamma squared.

    Given: two el. with gamma differing by a factor of 2
    When: synchrotron loss rate is computed for both
    Then: loss rate ratio must equal the square of the gamma ratio
    """
    loss_1 = synchrotron_loss(gamma=1e4, B=1e-9)
    loss_2 = synchrotron_loss(gamma=2e4, B=1e-9)
    assert np.isclose(loss_2 / loss_1, 4.0, rtol=1e-5)


def test_synchrotron_loss_scales_with_B_squared():
    """Test that synchrotron loss scales as B squared via U_B = B^2 / 2mu_0.

    Given: two identical el. in magnetic fields differing by a factor of 2
    When: synchrotron loss rate is computed for both
    Then: loss rate ratio must equal the square of the field ratio
    """
    loss_1 = synchrotron_loss(gamma=1e4, B=1e-9)
    loss_2 = synchrotron_loss(gamma=1e4, B=2e-9)
    assert np.isclose(loss_2 / loss_1, 4.0, rtol=1e-5)


def test_synchrotron_loss_is_zero_for_zero_field():
    """Test that synchrotron loss vanishes When there is no magnetic field.

    Given: an el. with a typical Lorentz factor and zero magnetic field
    When: the synchrotron loss rate is computed
    Then: result must be zero since B=0
    """
    assert synchrotron_loss(gamma=1e4, B=0.0) == 0.0


"///////////////////////////////////////////////////////////////////////////////////"
"///////////////////////////////////////////////////////////////////////////////////"


"COMPTON LOSS TESTS"
 


def test_inverse_compton_loss_is_negative():
    """Test that inverse Compton loss always removes energy from the el..

    Given: an el. with a typical Lorentz factor at redshift zero
    When: the inverse Compton loss rate is computed
    Then: the result must be negative since el. only lose energy
    """
    assert inverse_compton_loss(gamma=1e4, redshift=0.0) < 0


def test_inverse_compton_loss_scales_with_gamma_squared():
    """Test that inverse Compton loss scales as gamma squared.

    Given: two el. with Lorentz factors differing by a factor of 2
    When: inverse Compton loss rate is computed for both
    Then: loss rate ratio must equal the square of the gamma ratio
    """
    loss_1 = inverse_compton_loss(gamma=1e4, redshift=0.0)
    loss_2 = inverse_compton_loss(gamma=2e4, redshift=0.0)
    assert np.isclose(loss_2 / loss_1, 4.0, rtol=1e-5)


def test_inverse_compton_loss_increases_with_redshift():
    """Test that inverse Compton loss increases with redshift.

    Given: two identical el. at different redshifts
    When: inverse Compton loss rate is computed for both
    Then: loss at higher redshift must be larger because U_rad scales as (1+z)^4
    """
    loss_z0 = inverse_compton_loss(gamma=1e4, redshift=0.0)
    loss_z1 = inverse_compton_loss(gamma=1e4, redshift=1.0)
    assert abs(loss_z1) > abs(loss_z0)


def test_inverse_compton_loss_scales_with_redshift():
    """Test that inverse Compton loss scales as (1+z)^4 with redshift.

    Given: two identical el. at redshift 0 and redshift 1
    When: inverse Compton loss rate is computed for both
    Then: loss rate ratio must equal (1+1)^4 / (1+0)^4 = 16
    """
    loss_z0 = inverse_compton_loss(gamma=1e4, redshift=0.0)
    loss_z1 = inverse_compton_loss(gamma=1e4, redshift=1.0)
    assert np.isclose(loss_z1 / loss_z0, 16.0, rtol=1e-5)


"///////////////////////////////////////////////////////////////////////////////////"
"///////////////////////////////////////////////////////////////////////////////////"


"BREMMSTRAHLUNG LOSS TEST"

def test_bremsstrahlung_loss_is_negative():
    """Test that bremsstrahlung loss always removes energy from the el..

    Given: an el. with a typical Lorentz factor in a typical plasma
    When:  bremsstrahlung loss rate is computed
    Then:  result must be negative since el. only lose energy
    """
    assert bremsstrahlung_loss(gamma=1e4, n_plasma=1e3) < 0


def test_bremsstrahlung_loss_scales_with_gamma():
    """Test that bremsstrahlung loss scales linearly with gamma.

    Given: two el. with Lorentz factors differing by a factor of 2
    When:  bremsstrahlung loss rate is computed for both
    Then:  loss rate ratio must equal the gamma ratio
    """
    loss_1 = bremsstrahlung_loss(gamma=1e4, n_plasma=1e3)
    loss_2 = bremsstrahlung_loss(gamma=2e4, n_plasma=1e3)
    assert np.isclose(loss_2 / loss_1, 2.0, rtol=1e-5)


def test_bremsstrahlung_loss_scales_with_n_plasma():
    """Test that bremsstrahlung loss scales linearly with plasma density.

   Given: two identical el. in plasmas with densities differing by factor 2
 When:  bremsstrahlung loss rate is computed for both
    Then:  loss rate ratio must equal the density ratio
    """
    loss_1 = bremsstrahlung_loss(gamma=1e4, n_plasma=1e3)
    loss_2 = bremsstrahlung_loss(gamma=1e4, n_plasma=2e3)
    assert np.isclose(loss_2 / loss_1, 2.0, rtol=1e-5)


def test_bremsstrahlung_loss_is_zero_for_zero_density():
    """Test that bremsstrahlung loss vanishes When there is no plasma.

   Given: an el. with a typical Lorentz factor in an empty plasma
 When:  bremsstrahlung loss rate is computed
    Then:  result must be zero since there are no ions to scatter off
    """
    assert bremsstrahlung_loss(gamma=1e4, n_plasma=0.0) == 0.0



"///////////////////////////////////////////////////////////////////////////////////"
"///////////////////////////////////////////////////////////////////////////////////"


"COULOMB LOSS TEST"


def test_coulomb_loss_is_negative():
    """Test that Coulomb loss always removes energy from the el..

    Given: an el. with a typical Lorentz factor in a typical plasma
    When:  Coulomb loss rate is computed
    Then:  result must be negative since el. only lose energy
    """
    assert coulomb_loss(gamma=1e4, n_plasma=1e3) < 0


def test_coulomb_loss_scales_inversely_with_gamma():
    """Test that Coulomb loss scales as 1/gamma.

    Given: two el. with Lorentz factors differing by a factor of 2
    When:  Coulomb loss rate is computed for both
    Then:  loss rate ratio must equal 0.5 since loss goes as 1/gamma
    """
    loss_1 = coulomb_loss(gamma=1e4, n_plasma=1e3)
    loss_2 = coulomb_loss(gamma=2e4, n_plasma=1e3)
    assert np.isclose(loss_2 / loss_1, 0.5, rtol=1e-5)


def test_coulomb_loss_scales_with_n_plasma():
    """Test that Coulomb loss scales linearly with plasma density.

    Given: two identical el. in plasmas with densities differing by factor 2
    When:  Coulomb loss rate is computed for both
    Then:  loss rate ratio must equal the density ratio
    """
    loss_1 = coulomb_loss(gamma=1e4, n_plasma=1e3)
    loss_2 = coulomb_loss(gamma=1e4, n_plasma=2e3)
    assert np.isclose(loss_2 / loss_1, 2.0, rtol=1e-5)


def test_coulomb_loss_is_zero_for_zero_density():
    """Test that Coulomb loss vanishes When there is no plasma.

    Given: an el. with a typical Lorentz factor in an empty plasma
    When:  Coulomb loss rate is computed
    THEN:  result must be zero since there are no particles to collide with
    """
    assert coulomb_loss(gamma=1e4, n_plasma=0.0) == 0.0










