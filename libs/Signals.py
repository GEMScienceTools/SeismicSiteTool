# -*- coding: utf-8 -*-
#
# Copyright (C) 2010-2016 GEM Foundation
#
# The Site Response Toolkit (SRTK) is free software: you can redistribute
# it and/or modify it under the terms of the GNU Affero General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# SRTK is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# with this download. If not, see <http://www.gnu.org/licenses/>
#
# Author: Poggi Valerio

#-----------------------------------------------------------------------------------------

class Signal():

  def __init__(self):
    """
    HDR (Header, floats):
      NSMP - Number of samples
      TSMP - Sampling time
      TREF - Reference time: [Year, Month, Day, Hour, Minute, Second]
      SCRD - Station coordinates: [Lon., Lat., Depth] or [X, Y, Z]
      NCHN - Number of channels
      ROTA - Rotation angles of each channel: [α, β, γ, ...]

    INF (Meta Information, strings):
      NETID - Network ID
      STAID - Station ID
      RECID - Recording ID
      EVEID - Event ID
      CHNID - Channel ID (one name for each channel)
      UNITS - Units in S.I. or description (e.g. 'Acceleration')

    CHN - Recording list (one array per channel), floats
    TAX - Time axis of the recording (Optional), floats
    FAX - Frequency axis of the spectrum (Optional), floats
    """

    self.HDR = {'NSMP': [],
                'TSMP': [],
                'TREF': [],
                'SCRD': [],
                'NCHN': [],
                'ROTA': []}

    self.INF = {'NETID': '',
                'STAID': '',
                'RECID': '',
                'EVEID': '',
                'CHNID': '',
                'UNITS': ''}

    self.TAX = []
    self.FAX = []
    self.CHN = []

