#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright (C) 2010-2017 GEM Foundation
#
# The Site Response Toolkit (SRTK) is free software: you can redistribute
# it and/or modify it under the terms of the GNU Affero General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# SRTK is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# Author: Poggi Valerio

"""

"""

import numpy as _np

#-----------------------------------------------------------------------------------------

def IsEmpty(Number):
  """
  Check if variable is empty with different formats.
  """ 

  C0 = (Number == [])
  C1 = (Number == '')
  C2 = (Number != Number)
  C3 = (Number == None)
  C4 = (Number == 'None')

  Eval = (C0 or C1 or C2 or C3 or C4)

  return Eval

#-----------------------------------------------------------------------------------------

def NoneCheck(Number):
  """
  """

  Eval = IsEmpty(Number)
  Number = None if Eval else Number

  return Number

#-----------------------------------------------------------------------------------------

def Round(Number, Decimal=3):
  """
  Round scalar and arrays to a given decimal place
  """

  if isinstance(Number, (list, tuple, _np.ndarray)):
    for I,N in enumerate(Number):
      Number[I] = round(N, Decimal)
  else:
    Number = round(Number, Decimal)

  return Number
