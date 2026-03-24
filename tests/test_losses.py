"""
Tests for energy loss functions in plasma_sim.losses.

Each test verifies a single physical property of the loss functions.
"""

import numpy as np
from plasma_sim.losses import synchrotron_loss


def test_synchrotron_loss_is_negative():
    """Test that synchrotron loss always removes energy from the electron.

    Given: an electron with a typical Lorentz factor and a typical magnetic field
    When: the synchrotron loss rate is computed
    Then: the result must be negative since electrons only lose energy
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