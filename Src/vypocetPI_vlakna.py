import sys

import random
from PyQt5.QtCore import QThread, pyqtSignal, pyqtSlot
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QPushButton, QLabel, QWidget, QHBoxLayout, QMainWindow


class ComputationalPI(QThread):
    """
    Třída pro výpočet čísla pi pomocí Monte Carlo metody.
    Funguje jako služba pro výpočet, která je řízena GUI
    """
    progress_updated = pyqtSignal(int)  # Signal pro aktualizaci GUI
    value_updated = pyqtSignal(float)  # Signal pro aktualizaci GUI
    finished = pyqtSignal()  # Signal pro dokončení výpočtu

    def __init__(self, steps=100, one_step_iterations=10000):
        super().__init__()
        self.steps = steps  # Počet kroků výpočtu
        self.one_step_iteration = one_step_iterations  # Počet iterací v jednom kroku
        self.current_step = 0 # Aktuální krok výpočtu
        # atributy pro ovládání vlákna
        self._is_paused = False
        self._is_stopped = False
        # zde bude uložena hodnota pi
        self.pi_value =0
        self.cnt_in = 0
        self.cnt_total = 0


    def compute_one_step(self):
        # Simulace jednoho kroku výpočtu
        for i in range(self.one_step_iteration):
            x = random.random()
            y = random.random()
            self.cnt_total+=1
            if x ** 2 + y ** 2 <= 1:
                self.cnt_in += 1
        self.pi_value = 4 * self.cnt_in / self.cnt_total
        QThread.msleep(50)  # Ješte to zpozdi, aby bylo videt, ze se něco děje

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
        self._is_stopped = True

    def run(self):
        """
        Spustí výpočet
        :return:
        """
        # nastavení stavovývh proměnných
        self._is_stopped = False
        self._is_paused = False
        # dokud neproběhne požadovaný počet kroků tak opakujeme výpočet
        while self.current_step < self.steps:
            if self._is_stopped:  # pokud bylo vlákno zastaveno, skonči
                self.finished.emit()
                return
            if self._is_paused: # pokud bylo vlákno pozastaveno, počkej
                QThread.msleep(100)  # Čekání, dokud není odblokováno
                continue
            # Provádíme krok výpočtu
            self.compute_one_step() # spočti jeden krok
            self.current_step += 1 # zvyš počet provedených kroků
            # řekni GUI, že se změnil stav
            self.progress_updated.emit(self.current_step)
            self.value_updated.emit(self.pi_value)

        self.finished.emit() # po skončení výpočtu emituj signál, že je výpočet dokončen


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Výpočet \u03c0 ")
        self.setGeometry(100, 100, 300, 200)

        #============= Rozvržení aplikace ==============
        central_widget = QWidget()
        # Layout - vše umisťujeme ve sloupci
        vertical_layout = QVBoxLayout(central_widget)
        # Popisky - jsou umístěny vedle sebe v jednom řádku
        horizontal_layout = QHBoxLayout()
        self.progress_label = QLabel("Progres: 0%")
        horizontal_layout.addWidget(self.progress_label)
        self.result_label = QLabel("Pi: 0")
        horizontal_layout.addWidget(self.result_label)
        vertical_layout.addLayout(horizontal_layout)

        # --------------Přidání tlačítek----------------
        self.start_button = QPushButton("Start")
        self.start_button.setStyleSheet("background-color: green;font-weight: bold; font-size: 16px;")
        self.start_button.clicked.connect(self.start_computation)
        vertical_layout.addWidget(self.start_button)

        self.pause_button = QPushButton("Pause")
        self.pause_button.setStyleSheet("color: blue")
        self.pause_button.clicked.connect(self.pause_computation)
        vertical_layout.addWidget(self.pause_button)

        self.resume_button = QPushButton("Resume")
        self.resume_button.clicked.connect(self.resume_computation)
        vertical_layout.addWidget(self.resume_button)

        self.stop_button = QPushButton("Stop")
        self.stop_button.setStyleSheet("background-color: red")
        self.stop_button.clicked.connect(self.stop_computation)
        vertical_layout.addWidget(self.stop_button)
        # přidání vertikálního layoutu z central widgetu do hlavního okna
        self.setCentralWidget(central_widget)
        self.statusBar().showMessage("Připraveno")
        #===============================================

        self.worker = ComputationalPI(steps=100) # Vytvoření instance výpočetního vlákna
        self.worker.progress_updated.connect(self.update_progress_label) # Připojení slotu k signálu - update prograsu
        self.worker.value_updated.connect(self.update_pi_label)    # Připojení slotu k signálu - update hodnoty
        self.worker.finished.connect(self.computation_finished) # Připojení slotu k signálu - konec výpočtu

    #------------------------ Sloty pro ovládání výpočtu---------------------
    def start_computation(self):
        """
        Spustí výpočet
        :return:
        """
        self.statusBar().showMessage("Výpočet probíhá...")
        self.worker.current_step = 0 # vynuluje stav vlákna
        if not self.worker.isRunning():  # Spustíme vlákno pouze, pokud neběží
            self.worker.start()

    def pause_computation(self):
        """
        Pozastaví výpočet
        :return:
        """
        self.worker.pause()
        self.statusBar().showMessage("Výpočet pozastaven...")

    def resume_computation(self):
        """
        Odpauzuje výpočet
        :return:
        """
        self.worker.resume()
        self.statusBar().showMessage("Výpočet obnoven...")

    def stop_computation(self):
        self.worker.stop()
        self.statusBar().showMessage("Výpočet přerušen...")

    #------------------------ Sloty pro aktualizaci GUI---------------------
    @pyqtSlot(int)
    def update_progress_label(self, step):
        """
        Metoda pro aktualizaci popisku progresu
        :param step:
        :return:
        """
        progress = int((step / self.worker.steps) * 100)
        self.progress_label.setText(f"Progres: {progress}%")

    @pyqtSlot(float)
    def update_pi_label(self, value):
        """
        Metoda pro aktualizaci popisku s výsledkem
        :param value: 
        :return: 
        """
        self.result_label.setText(f"Pi: {value:.6f}")

    @pyqtSlot()
    def computation_finished(self):
        """
        Metoda, která se zavolá po dokončení výpočtu
        :return: 
        """
        self.statusBar().showMessage("Výpočet dokončen!")


if __name__ == "__main__":
    """
    Hlavní funkce programu
    Vytvoření instance QApplication, instance hlavního okna a spuštění aplikace
    """
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
    
