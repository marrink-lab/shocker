#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  9 11:53:09 2022

@author: marco
"""

try:
    import pkg_resources
except ImportError:
    import os
    DATA_PATH = os.path.join(os.path.dirname(__file__), 'data')
    TEST_DATA = os.path.join(os.path.dirname(__file__), 'tests/test_data')
    del os
else:
    DATA_PATH = pkg_resources.resource_filename('shocker', 'data')
    TEST_DATA = pkg_resources.resource_filename('shocker', 'tests/test_data')
    del pkg_resources