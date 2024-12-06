import sys

import matplotlib.pyplot as plt
import numpy as np
from PyQt5.QtCore import QThread, pyqtSignal, pyqtSlot
from PyQt5.QtGui import QImage, QPixmap

from LBM import LBM # import třídy pro výpočet LBM


class ComputationalLBM(QThread):
    """
    Třída pro simulaci proudění
    Funguje jako adaptér mezi GUI a LBM
    """
    progress_updated = pyqtSignal(float)  # Signal pro aktualizaci GUI - aktuální čas
    velocity_field_updated = pyqtSignal(QImage)  # Signal pro aktualizaci GUI - pole rychlostí
    finished = pyqtSignal()  # Signal pro dokončení výpočtu

    def __init__(self, T, u0,nu,R):
        super().__init__()
        # nastavení některých parametrů simulace na pevno, neumožňuji měnit
        Lx = 2.2 # rozměry oblasti - délka
        Ly = 0.41 # rozměry oblasti - šířka
        dx = 0.01 # velikost buňky
        self.T = T  # maximální doba výpočtu - simulovaný čas
        self.current_time = 0  # aktuální čas výpočtu
        # uMax, nu_phys, Lx, Ly, T, dx, Cx = 0.2, Cy = 0.2, R=0.05):
        self.lbm_instance = LBM(u0,nu, Lx, Ly,T, dx, R=R) # Vytvoření instance LBM

        # atributy pro ovládání vlákna
        self._is_paused = False
        self._is_stopped = False

    def pause(self):
        """
        Pozastaví výpočet
        :return:
        """
        self._is_paused = True

    def resume(self):
        """
        Odpauzuje výpočet
        :return:
        """
        self._is_paused = False

    def stop(self):
        """
        Zastaví výpočet
        :return:
        """
        self._is_stopped = True

    def run(self):
        """
        Spustí výpočet
        :return:
        """
        # nastavení stavových proměnných
        self._is_stopped = False
        self._is_paused = False
        self.lbm_instance.start()  # spust LBM
        while self.current_time<self.T: # dokud není dosaženo maximálního času
            if self._is_stopped: # pokud je zastaveno, skonči
                self.finished.emit() # emituj signál, že je výpočet dokončen
                return

            if self._is_paused: # pokud je pozastaveno, počkej
                QThread.msleep(100)  # Čekání, dokud není odblokováno
                continue
            # Provádíme krok výpočtu a emitujeme obrázek pro GUI
            self.lbm_instance.nextIteration() # proved krok simulace
            # po každém kroku simulace spočtu rychlost a udělám z toho obrázek
            # pozor - trochu prasárna, vytvářím spoustu instancí QImage, ale pro jednoduchost to nechám
            velocity = np.sqrt(self.lbm_instance.velocity2D[:,:,0]**2+self.lbm_instance.velocity2D[:,:,1]**2)
            normalized_velocity = ((velocity - np.min(velocity)) / (np.max(velocity) - np.min(velocity)) * 255).astype(np.uint8)
            image = QImage(normalized_velocity, velocity.shape[1], velocity.shape[0], QImage.Format_Grayscale8)
            self.current_time = self.lbm_instance.physicalTime()  # zjisti aktuální čas simulace
            # emituj signály pro aktualizaci GUI
            self.velocity_field_updated.emit(image) # emituj signál pro aktualizaci obrázku v GUI
            self.progress_updated.emit(self.current_time)

        # po skončení výpočtu emituj signál, že je výpočet dokončen
        self.finished.emit()

