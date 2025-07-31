import numpy as np
import pandas as pd
import csv
import time
import datetime
import serial
import struct
from pyModbusTCP.client import ModbusClient


class TemperatureModbus(object):
    """
    object of class TemperatureModbus interfaces CABTR (Centrale 
    d’Acquisition Basses Températures Rapide) - fast acquisition module
    for low temperatures and acquires up to 8 temperature values
    """
    def __init__(self, host='192.168.0.225', port=502, auto_open=True):
        try:
            self.c = ModbusClient()
            self.c.host(host)
            self.c.port(port)
            self.c.auto_open(True)
        except RuntimeError as e:
            raise RuntimeError('Runtime Error: {e}')
        except Exception as e:
            print("Unknown error: {e}")
            raise
    
    def __call__(self, i):
        """ return temperature in K for i-th channel"""
        assert i in range(1, 9)
        address = [272, 274, 276, 278, 280, 282, 284, 286][i-1]
        if self.c.read_holding_registers(166)[0]:
            f1 = self.c.read_holding_registers(address)[0]
            f2 = self.c.read_holding_registers(address+1)[0]
            T = struct.unpack('f', struct.pack('HH', *(f1, f2)))[0]
            return round(T, 8)
        else:
            raise ValueError('Temperature read error')


class PressureSerial(object):
    """
    object of class PressureSerial interfaces Mensor CPT 6100
    pressure sensors and allows to send queries for pressure
    values, as well as setting basic parameters
    """
    def __init__(self, port, address='*', timeout=2):
        self.address = address
        self.ser = serial.Serial()
        self.ser.port = port
        self.ser.timeout = timeout
        self.ser.baudrate = 9600 # bits/sec
        self.ser.bytesize = serial.EIGHTBITS
        self.ser.stopbits = serial.STOPBITS_ONE
        self.ser.parity = serial.PARITY_NONE
        try:
            self.ser.open()
        except:
            raise RuntimeError(f'{port} port error - unable to open')
    
    def setter(self, s):
        s = s + '\r'
        self.ser.write(s.encode())
    
    def query(self, s):
        s = s + '\r'
        self.ser.write(s.encode())
        return self.ser.readline().decode().split()
   
    def save(self, address=None):
        """
        save the current turn down data to EEPROM
        without executing, settings are valid during a single session
        """
        if address is None:
            address = self.address 
        self.setter('#{}SAVE'.format(address))

    def set_address(self, n, address=None):
        """
        permanently set sensor existing address to n (0-9 or A-Z)
        n is the new address
        address is the existing address
        """
        if address is None:
            address = self.address 
        self.setter('#{}A {}'.format(address, n))
        self.setter('#{}SAVE'.format(n))
        self.address = n
   
    def set_filter(self, value):
        """ set filter value ranging from 0 to 99 """
        assert value in range(0, 100)
        self.setter('#{}FL {}'.format(self.address, int(value)))
    
    def get_filter(self):
        """ return filter value ranging from 0 to 99 """
        return int(self.query('#{}FL?'.format(self.address))[2])
    
    def __call__(self):
        """ return pressure in MPa """
        res = float(self.query('#{}?'.format(self.address))[1]) # bar
        return round(res / 10, 8) # MPa


class GasAnalyzerSerial(object):
    """
    object of class GasAnalyzerSerial interfaces SRS BGA244HP gas analyzer
    and allows recording binary gas composition, as well as the speed of sound,
    uncertainty, and other parameters.
    """
    def __init__(self, port, timeout = 5):
        self.ser = serial.Serial()
        self.ser.port = port
        self.ser.timeout = timeout
        self.ser.baudrate = 9600 # bits/sec
        self.ser.bytesize = serial.EIGHTBITS
        self.ser.stopbits = serial.STOPBITS_ONE
        self.ser.parity = serial.PARITY_NONE
        self.ser.rtscts = True
        try:
            self.ser.open()
        except:
            raise RuntimeError('{} port error. Unable to open'.format(port))
        self.casDict = {'argon': '7440-37-1',
                        'helium': '7440-59-7',
                        'helium-4': '7440-59-7',
                        'hydrogen': '1333-74-0',
                        'methane': '74-82-8',
                        'nitrogen': '7727-37-9',
                        'neon': '7440-01-9',
                        'oxygen': '7782-44-7'}
    
    def setter(self, s):
        s = s + '\r\n'
        self.ser.write(s.encode())
    
    def query(self, s):
        s = s + '\r\n'
        self.ser.write(s.encode())
        return self.ser.readline().decode().strip('\r\n')
    
    def set_run_mode(self, i):
        """ set run mode to 1 for RUN or 0 for STOP """
        assert i in [0, 1]
        self.setter('RUNM{}'.format(i))
    
    def get_run_mode(self):
        """ get actual run mode (1 for RUN, 0 for STOP) """
        return int(self.query('RUNM?'))
    
    def set_primary_gas(self, fld):
        """ set primary gas """
        cas = self.casDict[fld.lower()]
        self.setter('GASP {}'.format(cas))
    
    def set_secondary_gas(self, fld):
        """ set secondary gas """
        cas = self.casDict[fld.lower()]
        self.setter('GASS {}'.format(cas))
    
    def swap_binary_gases(self):
        """ swap binary gases """
        self.setter('SWAP')
    
    def get_gauge_pressure(self):
        """ return pressure in MPa """
        return round(float(self.query('PRES? bar')) / 10, 10)
    
    def get_uncertainty(self):
        """
        return measurement uncertainty
        of the currently active instrument mode,
        valid for Binary Gas and Gas Purity Modes """
        return round(float(self.query('UNCT?%')) / 100, 10)
    
    def get_SOS(self):
        """ return speed of sound in m/s """
        return float(self.query('SSOS? m/s'))
    
    def get_cell_temp(self):
        """ return cell temperature in K """
        return float(self.query('TCEL? K'))
    
    def get_all(self):
        """
        return a list of ratio 1, ratio 2,
        gas temperature, analysis pressure,
        normalized speed of sound and block temperature
        in global units
        """
        return [float(i) for i in self.query('XALL?').split(',')]
   
    def __call__(self, i=1):
        """ return mole fraction of the i-th fluid, i in [1, 2] """
        assert i in [1, 2]
        return round(float(self.query('RATO? {}%'.format(i))) / 100, 10)


def acquireData(primary_gas,
                secondary_gas,
                t_acquisition_s,
                t_tempMemory_s,
                filename):
    """
    function to read all thermodynamic parameters of interest and save
    them to files in an infinite loop. Arguments:
        `primary_gas` is the first gas in binary mixture
        `secondary_gas` is the second gas in the binary mixture
        `t_acquisition_s` is the data acquisition interval in seconds
        `t_tempMemory_s` is the time interval for real-time plotting in seconds
    """
    assert t_acquisition_s >= t_tempMemory_s

    PT101 = PressureSerial('COM3', address='1')
    PT102 = PressureSerial('COM6', address='2')
    GA016 = GasAnalyzerSerial('COM5')
    MDB = TemperatureModbus()

    GA016.set_primary_gas(primary_gas)
    GA016.set_secondary_gas(secondary_gas)

    with open(filename, "w", newline="") as fh:
        fieldnames = ["date/yyyy-mm-dd",
                      "time/hh:mm:ss",
                      "PT102/MPa",
                      "PT101/MPa",
                      "TT009/K",
                      "TT010/K",
                      "TT008/K",
                      "TT101/K",
                      "TT102/K",
                      "TT006/K",
                      "TT007/K",
                      "x1",
                      "x2",
                      "err(x)",]
        csv_writer = csv.DictWriter(fh, fieldnames=fieldnames)
        csv_writer.writeheader()

        j = 0
        while True:
            # get date and time
            date, t = str(datetime.datetime.now()).split()

            # write data to dictionary
            data = {
                "date/yyyy-mm-dd": date,
                "time/hh:mm:ss": t.split(".")[0],
                "PT102/MPa": PT102(),
                "PT101/MPa": PT101(),
                "TT009/K": MDB(6),
                "TT010/K": MDB(8),
                "TT008/K": MDB(1),
                "TT101/K": MDB(5),
                "TT102/K": MDB(4),
                "x1": GA016(1),
                "x2": GA016(2),
                "err(x)": GA016.get_uncertainty(),
            }

            # dump to buffer
            with open("buffer.txt", "w", newline="") as textBuffer:
                textBuffer.truncate(0)
                pd.DataFrame(data, index=[0]).to_csv(textBuffer, index=False)

            # print header to console (TT006, TT007, TT008, TT009, TT010
            # are of secondary imporance and are not printed)
            if j % int((t_acquisition_s/t_tempMemory_s) * 10) == 0:
                print(["time/hh:mm:ss",
                       "PT101/MPa",
                       "PT102/MPa",
                       "TT101/K",
                       "TT102/K",
                       "x1",
                       "x2",
                       "err(x)",])

            if j % int(t_acquisition_s/t_tempMemory_s) == 0:
                # dump to file
                csv_writer.writerow(data)
                fh.flush()

                # print data to console
                print([data["time/hh:mm:ss"],
                       data["PT101/MPa"],
                       data["PT102/MPa"],
                       data["TT101/K"],
                       data["TT102/K"],
                       data["x1"],
                       data["x2"],
                       data["err(x)"]])
                
            j += 1
            time.sleep(t_tempMemory_s)
    PT01.close()
    DPT02.close()
    GA18.close()


def main():
    primary_gas = "nitrogen"
    secondary_gas = "helium"
    time_now = str(datetime.datetime.now().replace(microsecond=0)).replace(" ", "_")
    filename = f"data_{time_now}.csv"
    t_acquisition_s = 1 # sec
    t_tempMemory_s = 0.2 # sec

    acquireData(primary_gas = primary_gas,
                secondary_gas = secondary_gas,
                t_acquisition_s = t_acquisition_s,
                t_tempMemory_s = t_tempMemory_s,
                filename = filename)


if __name__ == "__main__":
    main()
