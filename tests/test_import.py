#!/usr/bin/env python

import pytest


def test_importability():
    from covid_serology import config

    assert config.paths is not None
