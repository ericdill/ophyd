# vi: ts=4 sw=4
'''
:mod:`ophyd.control.detector` - Ophyd Detectors Class
=====================================================

.. module:: ophyd.control.detector
   :synopsis:
'''

from __future__ import print_function
from .signal import (EpicsSignal, SignalGroup)


class DetectorStatus(object):
    def __init__(self, detector):
        self.done = False
        self.detector = detector
    def _finished(self, success=True, **kwargs):
        self.done = True



class Detector(SignalGroup):
    '''A Base Detector class

    Subclass from this to implement your own detectors
    '''

    #SUB_START = 'start_acquiring'
    SUB_DONE = 'done_acquiring'
    _SUB_REQ_DONE = '_req_done'  # requested move finished subscription
    def __init__(self, *args, **kwargs):
        super(Detector, self).__init__(*args, **kwargs)

    def configure(self, *args, **kwargs):
        '''Configure the detector for data collection.

        This method configures the Detector for data collection and is called
        before data collection starts.
        '''
        pass

    def deconfigure(self):
        '''Unset configuration of Detector

        This method resets the Detector and is called after data collection
        has stopped.
        '''
        pass

    def acquire(self, **kwargs):
        '''Start an acquisition on the detector (c.f. Trigger)

        This routine starts a data acquisition and returns an object which is
        the status of the acquisition.

        Returns
        -------
        DetectorStatus : Object to tell if detector has finished acquiring
        '''
        status = DetectorStatus(self)
        self.subscribe(status._finished,
                       event_type=self._SUB_REQ_DONE, run=False)
        return status

    def _done_acquiring(self, timestamp=None, value=None, **kwargs):
        '''Call when acquisition has completed.  Runs SUB_DONE subscription.'''

        self._run_subs(sub_type=self.SUB_DONE, timestamp=timestamp,
                       value=value, **kwargs)

        self._run_subs(sub_type=self._SUB_REQ_DONE, timestamp=timestamp,
                       value=value, success=True,
                       **kwargs)
        self._reset_sub(self._SUB_REQ_DONE)

    def read(self, **kwargs):
        '''Retrieve data from instrumentation, format it, and return it.
        '''
        raise NotImplementedError('Detector.read must be implemented')


class SignalDetector(Detector):
    def __init__(self, signal, *args, **kwargs):
        super(SignalDetector, self).__init__(*args, **kwargs)
        self._signal = signal

    def read(self, **kwargs):
        '''Read signal and return formatted for run-engine.

        Returns
        -------
        dict
        '''
        return {self._signal.name: {'value': self._signal.value,
                                    'timestamp': self._signal.timestamp}}