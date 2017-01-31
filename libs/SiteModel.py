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
Base class to store and a analyse site information
"""

import numpy as _np
import SiteMethods as _SM
import AsciiTools as _AT

#-----------------------------------------------------------------------------------------

HdrKeys = ['Id','X','Y','Z']
ModKeys = ['Hl','Vp','Vs','Dn','Qp','Qs']
EngKeys = ['Vz','Qwl','K0','EC8']
AmpKeys = ['Freq','Shtf','Qwl','Attf','Fn','An']

#-----------------------------------------------------------------------------------------

def _IsEmpty(Number):
  """
  INTERNAL: Checking if variable is empty. Return boolean [0] and None [1]
  """ 

  C0 = (Number == [])
  C1 = (Number == '')
  C2 = (Number != Number)
  C3 = (Number == None)
  C4 = (Number == 'None')

  Eval = (C0 or C1 or C2 or C3 or C4)
  Number = None if Eval else Number

  return (Eval, Number)

#-----------------------------------------------------------------------------------------

class Site1D(object):
  """
  """

  def __init__(self, Id=[], X=[], Y=[], Z=[]):
    """
    """

    self.Hdr = {}
    self.Mod = {}
    self.Eng = {}
    self.Amp = {}

    for K in HdrKeys:
      self.Hdr[K] = []
    for K in ModKeys:
      self.Mod[K] = []
    for K in EngKeys:
      self.Eng[K] = []
    for K in AmpKeys:
      self.Amp[K] = []

    self.Hdr['Id'] = _IsEmpty(Id)[1]
    self.Hdr['X'] = _IsEmpty(X)[1]
    self.Hdr['Y'] = _IsEmpty(Y)[1]
    self.Hdr['Z'] = _IsEmpty(Z)[1]

  #---------------------------------------------------------------------------------------

  def AddLayer(self, Data, Index=-1):
    """
    Method to add a single layer to the site model, at the bottom of an existing stack
    (default) or at arbitrary location (Index). Data can be a list of values, sorted
    according to the ModKeys list, or a dictionary with corresponding keys.
    """

    Empty = None

    # Case: List
    if type(Data) is list:
      for I, K in enumerate(ModKeys):
        if Index == -1: Index = len(self.Mod[K])
        self.Mod[K].insert(Index, _IsEmpty(Data[I])[1])

    # Case: Dictionary
    if type(Data) is dict:
      for K in ModKeys:
        if Index == -1: Index = len(self.Mod[K])
        if K in Data.keys():
          self.Mod[K].insert(Index, _IsEmpty(Data[K])[1])

  #---------------------------------------------------------------------------------------

  def ImportModel(self, AsciiFile,
                        Header=[],
                        Delimiter=',',
                        Skipline=0,
                        Comment='#'):
    """
    Method to parse soil properties from an ascii file.
    Arbitray ascii format is allowed (default is CSV).
    """

    Table = _AT.AsciiTable()
    Table.Import(AsciiFile, header=Header,
                            dtype='float',
                            delimiter=Delimiter,
                            skipline=Skipline,
                            comment=Comment)

    for D in Table.data:
      self.AddLayer(D)

  #---------------------------------------------------------------------------------------

  def ComputeTTAV(self, Key='Vs', Z=30.):
    """
    Compute and store travel-time average velocity at a given depth (Z)
    and for a specific key. Default is Vs30. Multiple depths are also allowed.
    """

    if type(Z) != list:
      Z = [Z]

    # Initialise Vz data structure, if empty
    if _IsEmpty(self.Eng['Vz'])[0]:
      self.Eng['Vz'] = {}

    for z in Z:
      # Compute average velocity
      Vz = _SM.TTAverageVelocity(self.Mod['Hl'], self.Mod[Key], z)

      self.Eng['Vz']['{}'.format(z)] = Vz

  #---------------------------------------------------------------------------------------

  def ComputeQWL(self, Key='Vs'):
    """
    Compute and store the quarter-wavelength parameters.
    """

    # Initialise Qwl data structure, if empty
    if _IsEmpty(self.Eng['Qwl'])[0]:
      self.Eng['Qwl'] = {}

    # Compute average velocity
    Qwl = _SM.QwlApproxSolver(self.Mod['Hl'],
                              self.Mod[Key],
                              self.Mod['Dn'],
                              self.Amp['Freq'])

    self.Eng['Qwl']['Hl'] = Qwl[0]
    self.Eng['Qwl']['Vs'] = Qwl[1]
    self.Eng['Qwl']['Dn'] = Qwl[2]
    self.Amp['Qwl'] = Qwl[3]

  #---------------------------------------------------------------------------------------

  def ComputeKappa0(self, Key=('Vs','Qs'), Z=[]):
    """
    Compute the Kappa parameter from the Q profile of the site, down to
    a given depth (default is using the last interface of the profile).
    Multiple depths are not supported.
    """

    # Compute attenuation
    K0 = _SM.Kappa0(self.Mod['Hl'],
                    self.Mod[Key[0]],
                    self.Mod[Key[1]],
                    Z)

    self.Eng['K0'] = K0

  #---------------------------------------------------------------------------------------

  def ComputeGTClass(self, BCode='EC8'):
    """
    Compute geotechnical classification according to specified building code.
    Default is EC8. (Missing special classes)
    """

    try:
      Vs30 = self.Eng['Vz']['30.0']
    except:
      print 'Warning: Vs30 not found'
      return

    # Check the values
    if Vs30 >= 800.:
      self.Eng['EC8'] = 'A'
    if Vs30 >= 360. and Vs30 < 800.:
      self.Eng['EC8'] = 'B'
    if Vs30 >= 180. and Vs30 < 360.:
      self.Eng['EC8'] = 'C'
    if Vs30 < 180.:
      self.Eng['EC8'] = 'D'

  #---------------------------------------------------------------------------------------

  def ComputeSHTF(self, Iang=0., Elastic=False):
    """
    Compute the SH transfer function for an arbitrary incidence angle.
    Default incidence is vertical.
    """

    # TF calculation
    Shtf = _SM.ShTransferFunction(self.Mod['Hl'],
                                  self.Mod['Vs'],
                                  self.Mod['Dn'],
                                  self.Mod['Qs'],
                                  self.Amp['Freq'],
                                  Iang, Elastic)

    self.Amp['Shtf'] = Shtf

  #---------------------------------------------------------------------------------------

  def ComputeFnRes(self):
    """
    Identify resonance frequencies of the SH-wave transfer function.
    """

    # Get frequencies
    Fn, An = _SM.GetResFreq(self.Amp['Freq'], self.Amp['Shtf'])

    self.Amp['Fn'] = Fn
    self.Amp['An'] = An

  #---------------------------------------------------------------------------------------

  def ComputeAttFun(self):
    """
    Compute the frequency-dependent attenuation function for a given Kappa0.
    """

    # Compute exponential decay function
    Attf = _SM.AttenuationDecay(self.Amp['Freq'], self.Eng['K0'])

    self.Amp['Attf'] = Attf

  #---------------------------------------------------------------------------------------

  def GenFreqAx(self, Fmin=0.1, Fmax=10., Fnum=100, Log=True):
    """
    Method to generate a lin/log spaced frequency axis.
    """

    # Generate frequency axis
    Freq = _SM.FrequencyAxis(Fmin, Fmax, Fnum, Log)

    self.Amp['Freq'] = Freq

#-----------------------------------------------------------------------------------------

class SiteDb(object):
  """
  """

  def __init__(self, Id=[], Name=[]):
    """
    """

    self.Info = {'Id': Id, 'Name': Name}
    self.Site = []

  #---------------------------------------------------------------------------------------

  def AddSite(self, Site=[]):
    """
    """
    if Site:
      self.Site.append(Site)
    else:
      self.Site.append(Site1D())

  #---------------------------------------------------------------------------------------

  def Size(self):
    """
    Method to return size of the database.
    """

    return len(self.Sites)

