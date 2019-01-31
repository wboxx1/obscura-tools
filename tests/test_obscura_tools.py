#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `obscura_tools` package."""

import pytest
import numpy as np

import obscura_tools.radome_obscura as obscura


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_obscura_radome():
    ans = obscura.radome(60, -6, 20, 300)
    assert len(ans['azimuths']) == 39
