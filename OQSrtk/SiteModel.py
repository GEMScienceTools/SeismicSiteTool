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
Base class to store and analyse site information
"""

import numpy as _np

import SiteMethods as _SM
import AsciiTools as _AT
import Utils as _UT

#-----------------------------------------------------------------------------------------

# Precision for decimal rounding
Decimal = 4

#-----------------------------------------------------------------------------------------

class Model(object):
  """
  Base class to store a single site model, including the vertical
  soil profile and derived engineering parameters.
  """

  ParKeys = ['Hl','Vp','Vs','Dn','Qp','Qs']
  EngKeys = ['Vz','Qwl','K0','Gc']
  AmpKeys = ['Stf','Imp','Att','Res']

  #---------------------------------------------------------------------------------------

  def __init__(self):

    self.ParInit()
    self.EngInit()
    self.AmpInit()

  def ParInit(self):
    self.Par = {}
    for K in self.ParKeys: self.Par[K] = []

  def EngInit(self):
    self.Eng = {}
    for K in self.EngKeys: self.Eng[K] = []

  def AmpInit(self):
    self.Amp = {}
    for K in self.AmpKeys: self.Amp[K] = []

  #---------------------------------------------------------------------------------------

  def AddLayer(self, Data, Index=-1):
    """
    Method to add a data layer to the soil profile at arbitrary location.
    Data can be a list of values, sorted according to the ModKeys list,
    or a dictionary with corresponding keys.
    """

    Index = int(Index)

    # Case: List
    if type(Data) is list:
      for I, K in enumerate(self.ParKeys):
        if Index < 0: Index = len(self.Par[K])
        self.Par[K].insert(Index, _UT.NoneCheck(Data[I]))

    # Case: Dictionary
    if type(Data) is dict:
      for K in self.ParKeys:
        if Index < 0: Index = len(self.Par[K])
        if K in Data.keys():
          self.Par[K].insert(Index, _UT.NoneCheck(Data[K]))

  #---------------------------------------------------------------------------------------

  def DelLayer(self, Index=-1):
    """
    """

    Index = int(Index)
    del self.Par[Index]

#-----------------------------------------------------------------------------------------

class Site1D(object):
  """
  """

  def __init__(self, Id=[], X=[], Y=[], Z=[]):

    self.Hdr = {}
    self.Hdr['Id'] = Id
    self.Hdr['X'] = X
    self.Hdr['Y'] = Y
    self.Hdr['Z'] = Z

    self.Mod = []

  #---------------------------------------------------------------------------------------

  def AddModel(self, Index=-1, Mod=[]):
    """
    Add a soil model to the site database
    """

    Index = int(Index)
    if Index < 0: Index = len(self.Mod)

    if Mod:
      self.Mod.insert(Index, Mod)
    else:
      self.Mod.insert(Index, Model())

  #---------------------------------------------------------------------------------------

  def DelModel(self, Index=-1):
    """
    Remove a soil model from the site database
    """

    Index = int(Index)
    del self.Mod[Index]

  #---------------------------------------------------------------------------------------

  def ImportModel(self, AsciiFile, FileType='',
                                   Header=[],
                                   Delimiter=',',
                                   SkipLine=0,
                                   Comment='#',
                                   Index=-1,
                                   Owrite=False):
    """
    Method to parse soil properties from a single ascii file.
    Arbitray ascii format is allowed (default is CSV with header).
    """

    if Owrite and not _UT.IsEmpty(self.Mod):
      self.Mod[Index].ParInit()
    else:
      self.AddModel(Index)

    # STANDARD: CSV without header
    if FileType == 'CSV-NoH':
      Header = self.Mod[Index].ParKeys
      Delimiter = ','
      SkipLine = 0

    # STANDARD: Standard Geopsy format
    if FileType == 'Geopsy':
      Header = self.Mod[Index].ParKeys
      Delimiter = ' '
      SkipLine = 1

    Table = _AT.AsciiTable()
    Table.Import(AsciiFile, header=Header,
                            dtype='float',
                            delimiter=Delimiter,
                            skipline=SkipLine,
                            comment=Comment)

    Index = int(Index)

    for D in Table.data:
      self.Mod[Index].AddLayer(D)

  #---------------------------------------------------------------------------------------

  def FrequencyAxis(self, Fmin=0.1, Fmax=100., Fnum=1000, Log=True):
    """
    Method to generate a lin/log spaced frequency axis
    """
    self.Freq = _SM.FrequencyAxis(Fmin, Fmax, Fnum, Log)

  #---------------------------------------------------------------------------------------

  def ComputeTTAV(self, Key='Vs', Z=30.):
    """
    Compute and store travel-time average velocity at a given depth (Z).
    Default is Vs30. Multiple depths are also allowed (as list).
    """

    if type(Z) != list:
      Z = [Z]

    for M in self.Mod:

      # Initialise Vz data structure
      M.Eng['Vz'] = {}

      for z in Z:
        # Compute average velocity
        Vz = _SM.TTAverageVelocity(M.Par['Hl'], M.Par[Key], z)

        M.Eng['Vz'][z] = _UT.Round(Vz, Decimal)

  #---------------------------------------------------------------------------------------

  def ComputeGTClass(self, BCode='EC8'):
    """
    Compute geotechnical classification according to specified building code.
    Default is EC8. (Missing special classes)
    """

    for M in self.Mod:

      try:
        Vs30 = M.Eng['Vz'][30.]
      except:
        self.ComputeTTAV()
        Vs30 = M.Eng['Vz'][30.]

      if BCode == 'EC8':
        if Vs30 >= 800.:
          M.Eng['Gc'] = 'A'
        if Vs30 >= 360. and Vs30 < 800.:
          M.Eng['Gc'] = 'B'
        if Vs30 >= 180. and Vs30 < 360.:
          M.Eng['Gc'] = 'C'
        if Vs30 < 180.:
          M.Eng['Gc'] = 'D'

  #---------------------------------------------------------------------------------------

  def ComputeQWL(self, Key='Vs'):
    """
    Compute quarter-wavelength parameters and store
    them into the site database
    """

    for M in self.Mod:

      # Compute average velocity
      Qwl = _SM.QwlApproxSolver(M.Par['Hl'],
                                M.Par[Key],
                                M.Par['Dn'],
                                self.Freq)

      M.Eng['Qwl'] = {}
      M.Eng['Qwl']['Hl'] = _UT.Round(Qwl[0], Decimal)
      M.Eng['Qwl'][Key] = _UT.Round(Qwl[1], Decimal)
      M.Eng['Qwl']['Dn'] = _UT.Round(Qwl[2], Decimal)

  #---------------------------------------------------------------------------------------

  def ComputeImpAmp(self, Key='Vs', Vref=[], Dref=[]):
    """
    Compute the impedance amplification from Qwl parameters
    """

    for M in self.Mod:

      if not Vref:
        Vref = M.Par[Key][-1]
      if not Dref:
        Dref = M.Par['Dn'][-1]

      Amp = _SM.QwlImpedance(M.Eng['Qwl'][Key],
                             M.Eng['Qwl']['Dn'],
                             Vref, Dref)

      M.Amp['Imp'] = _UT.Round(Amp, Decimal)

  #---------------------------------------------------------------------------------------

  def ComputeSHTF(self, Iang=0., Elastic=False):
    """
    Compute the SH transfer function for an arbitrary incidence angle.
    Default incidence is vertical.
    """

    for M in self.Mod:

      # TF calculation
      Shtf = _SM.ShTransferFunction(M.Par['Hl'],
                                    M.Par['Vs'],
                                    M.Par['Dn'],
                                    M.Par['Qs'],
                                    self.Freq,
                                    Iang, Elastic)

      M.Amp['Stf'] = Shtf

  #---------------------------------------------------------------------------------------

  def ComputeFnRes(self):
    """
    Identify resonance frequencies of the SH-wave transfer function.
    """

    for M in self.Mod:

      try:
        Shtf = M.Amp['Stf']
      except:
        print 'Warning: Transfer Function not found'
        return

      # Get frequencies
      Fn, An = _SM.GetResFreq(self.Freq, Shtf)

      M.Amp['Res'] = {}
      M.Amp['Res']['Fn'] = _UT.Round(Fn, Decimal)
      M.Amp['Res']['An'] = _UT.Round(An, Decimal)

  #---------------------------------------------------------------------------------------

  def ComputeKappa(self, Key=('Vs','Qs'), Z=[]):
    """
    Compute the Kappa parameter from the Q profile of the site, down to
    a given depth (default is using the last interface of the profile).
    Multiple depths are not supported.
    """

    for M in self.Mod:

      # Compute kappa attenuation
      K0 = _SM.Kappa0(M.Par['Hl'],
                      M.Par[Key[0]],
                      M.Par[Key[1]],
                      Z)

      M.Eng['K0'] = _UT.Round(K0, Decimal)

  #---------------------------------------------------------------------------------------

  def ComputeAttFun(self):
    """
    Compute the frequency-dependent attenuation function for a given Kappa0.
    """

    for M in self.Mod:

      try:
        Kappa0 = M.Eng['K0']
      except:
        print 'Warning: Kappa 0 not found'
        return

      # Compute exponential decay function
      Attf = _SM.AttenuationDecay(self.Freq, Kappa0)

      M.Amp['Att'] = _UT.Round(Attf, Decimal)

#-----------------------------------------------------------------------------------------

class SiteDb(object):
  """
  """

  def __init__(self, Id=[], Info=[]):
    """
    """

    self.Hdr = {}
    self.Hdr['Id'] = Id
    self.Hdr['Info'] = Info

    self.Site = []

  #---------------------------------------------------------------------------------------

  def AddSite(self, Index=-1, Site=[]):
    """
    """

    Index = int(Index)
    if Index < 0: Index = len(self.Site)

    if Site:
      self.Site.insert(Index, Site)
    else:
      self.Site.insert(Index, Site1D())

  #---------------------------------------------------------------------------------------

  def ImportSites(self, AsciiFile, Root='', Filetype=''):
    """
    Import multiple site models from a csv list with format:
    'Id','X','Y','Z','File'
    """

    Table = _AT.AsciiTable()
    Table.Import(AsciiFile, header=['Id','X','Y','Z','File'],
                            dtype=['string','float','float','float','string'],
                            delimiter=',', 
                            skipline=0,
                            comment='#')

    for D in Table.data:
      S = Site1D(Table.data['Id'],
                 Table.data['X'],
                 Table.data['Y'],
                 Table.data['Z'])
      S.ImportModel(Table.data['File'],
                    Filetype=Filetype)
      self.AddSite(S)

  #---------------------------------------------------------------------------------------

  def ComputeTTAV(self, Key='Vs', Z=30., Average=False):
    """
    Compute average velocities for all sites in the database.
    """

    if type(Z) != list:
      Z = [Z]

    """
    # TO CHECK
    # Initialise Vz data structure, if empty
    if _IsEmpty(self.Eng['Vz'])[0]:
      self.Eng['Vz'] = {}

    for S in self.Site:
      S.ComputeTTAV(Key=Key, Z=Z)

    if Average:
      for z in Z:
        Data = [S.Eng['Vz'][z] for S in self.Site]
        Mn = _np.exp(_np.mean(_np.log(Data)))
        Sd = _np.exp(_np.std(_np.log(Data)))
        self.Eng['Vz'][float(z)] = [Mn, Sd]
    """

  #---------------------------------------------------------------------------------------

  def Size(self):
    """
    Method to return size of the database.
    """

    return len(self.Site)

