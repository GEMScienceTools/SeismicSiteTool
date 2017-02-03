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

import numpy as _np
import scipy.signal as _sig

import SacLib as _SL

#-----------------------------------------------------------------------------------------

class Record():

  def __init__(self):
    """
    HDR (Header, floats):
      NSMP - Number of samples
      TSMP - Sampling time
      TREF - Reference time: [Year, Month, Day, Hour, Minute, Second]
      SCRD - Station coordinates: [Lon., Lat., Depth/Elevation] or [X, Y, Z]
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

    self.HDR = {'NSMP': 0.,
                'TSMP': 0.,
                'TREF': [0,0,0,0,0,0.],
                'SCRD': [0,0,0],
                'NCHN': 0,
                'ROTA': [0,0,0]}

    self.INF = {'NETID': '',
                'STAID': '',
                'RECID': '',
                'EVEID': '',
                'CHNID': '',
                'UNITS': ''}

    self.TAX = []
    self.FAX = []
    self.CHN = []
    self.FSP = []

  #---------------------------------------------------------------------------------------

  def ImportSac(self, SacFile):
    """
    """

    S = _SL.Sac(SacFile)

    self.HDR['NSMP'] = S.Head['NPTS']
    self.HDR['TSMP'] = S.Head['DELTA']
    self.HDR['NCHN'] += 1

    self.CHN.append(_np.array(S.Data[0]))

  #---------------------------------------------------------------------------------------

  def Fourier(self, Inverse=False):
    """
    FFT
    """

    if Inverse:
      self.CHN = []
      for S in self.FSP:
        self.CHN.append(_np.fft.ifft(S).real)
    else:
      self.FSP = []
      for S in self.CHN:
        self.FSP.append(_np.fft.fft(S))

  #---------------------------------------------------------------------------------------

  def Taper(self, Alpha):
    """
    Tukey window tapering
    """

    Win = _sig.tukey(self.HDR['NSMP'], alpha=Alpha)

    for I,S in enumerate(self.CHN):
      self.CHN[I] *= Win

  #---------------------------------------------------------------------------------------

  def Filter(self, LowCorner, HighCorner, Order=3):
    """
    Butterworth bandpass filter
    """

    FS = 1./self.HDR['TSMP']

    if HighCorner >= FS/2.:
      print 'Warning: High corner must be < {0:.2f} Hz'.format(FS/2.)
      return

    if LowCorner < 0.:
      print 'Warning: Low corner must be > 0 Hz'.format(FS/2.)
      return

    # Corner frequencies
    Corners = [2.*LowCorner/FS, 2.*HighCorner/FS]

    # Butterworth filter
    b, a = _sig.butter(Order, Corners, btype='band')

    # Filtering records
    for I,S in enumerate(self.CHN):
      # self.CHN[I] = _sig.lfilter(b, a, S)
      zi = _sig.lfilter_zi(b, a);
      self.CHN[I],_ = _sig.lfilter(b, a, S, zi=zi*S[0])

#-----------------------------------------------------------------------------------------

class Array(object):
  """
  """

  def __init__(self, Id=[], Name=[]):
    """
    """

    self.Info = {'Id': Id, 'Name': Name}
    self.Record = []

  #---------------------------------------------------------------------------------------

  def AddRecord(self, Station=[]):
    """
    """
    if Station:
      self.Record.append(Station)
    else:
      self.Record.append(Record())
