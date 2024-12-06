import sys

from PyQt5.QtCore import Qt
from PyQt5.QtCore import pyqtSlot
from PyQt5.QtGui import  QImage, QPixmap
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QPushButton, QLabel, QWidget, QHBoxLayout, QMainWindow, \
    QFormLayout, QSpinBox, QDoubleSpinBox, QSlider

from Src.computationalLBM import ComputationalLBM


class SimulationGui(QMainWindow):
    """
    Třída pro simulaci proudění
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simulace proudění")
        self.setGeometry(100, 100, 840, 480) # Nastavení velikosti okna
        self.ratio = 3 # poměr jak se má přepočítat velikost obrázku na pixely
        self.image_size = (int(220*self.ratio), int(41*self.ratio))
        self.setup_gui()
        # nastaveni defaultnich hodnot pro simulace, dojde k jejich změně při změně hodnot v GUI
        self.T = 100
        self.R = 0.1
        self.nu = 0.1
        self.V0 = 1
        # nastavi defaltni hodnoty parametru, dle konfigurace GUI
        self.update_T()
        self.update_nu()
        self.update_V0()
        self.update_radius()

        # Přiřazení funkcí k tlačítkům
        self.start_button.clicked.connect(self.start_computation)
        self.pause_button.clicked.connect(self.pause_computation)
        self.resume_button.clicked.connect(self.resume_computation)
        self.stop_button.clicked.connect(self.stop_computation)
        # Přiřazení funkcí k signálům
        self.T_input.valueChanged.connect(self.update_T) # Přiřazení funkce k signálu změny hodnoty
        self.nu_input.valueChanged.connect(self.update_nu) # Přiřazení funkce k signálu změny hodnoty
        self.radius_input.valueChanged.connect(self.update_radius) # Přiřazení funkce k signálu změny hodnoty

        # Vytvoření workeru a přiřazení funkcí k ovládacím signálům
        #uMAx: 0.3, nu: 0.001, Lx: 2.2, Ly: 0.41, T: 1000, dx: 0.01
        self.worker = ComputationalLBM(T=self.T, u0=self.V0, nu=self.nu, R=self.R)  # Vytvoření instance výpočetního vlákna
        self.worker.progress_updated.connect(self.update_progress)  # Připojení slotu k signálu - update prograsu
        self.worker.velocity_field_updated.connect(self.update_image)  # Připojení slotu k signálu - update hodnoty
        self.worker.finished.connect(self.computation_finished)  # Připojení slotu k signálu - konec výpočtu
        # pomocná proměnná pro kontrolu, zda byly změny parametrů v GUI
        self.setting_updated = False # nastavení pro kontrolu, zda byly změny parametrů v GUI

    def setup_gui(self):
        """
        Vytvoření GUI
        :return:
        """
        central_widget = QWidget()
        vertical_layout = QVBoxLayout(central_widget)
        # =======================================================================================================
        # vytvoření hlavního layout, který obsahuje dva podlayouty - levý a pravý panel
        horizontal_layout_main = QHBoxLayout()
        # ---------levý panel obsahující formulář-----------
        left_form_layout = QFormLayout()
        # přidání zadávacích polí
        #  zadání doby běhu simulace
        self.T_label = QLabel("T - celkový čas")
        self.T_input = QSpinBox(self, minimum=1, maximum=10, value=1)
        left_form_layout.addRow(self.T_label, self.T_input)
        # zadání poloměru kruhu
        self.nu_label = QLabel("\u03BD - vyskozita")
        self.nu_input = QDoubleSpinBox(self, minimum=0.001, maximum=0.5, value=0.001, singleStep=0.0001)
        self.nu_input.setDecimals(4)
        left_form_layout.addRow(self.nu_label, self.nu_input)
        # zadání počáteční rychlosti
        self.V0_label = QLabel("V0 - počáteční rychlost")
        self.V0_input = QDoubleSpinBox(self, minimum=0.01, maximum=1, value=0.3, singleStep=0.1)
        left_form_layout.addRow(self.V0_label, self.V0_input)
        # zadání poloměru kruhu
        self.radius_label = QLabel("Poloměr:")
        self.radius_input = QSlider(Qt.Horizontal, self, minimum=1, maximum=10, value=5)
        left_form_layout.addRow(self.radius_label, self.radius_input)
        horizontal_layout_main.addLayout(left_form_layout)  # přidelíme levý panel do layoutu
        # ---------pravý panel obsahující obrázek-----------
        self.image_label = QLabel()
        self.image_label.setFixedSize(*self.image_size)
        self.image_label.setStyleSheet("border: 1px solid black")
        self.image_label.setAlignment(Qt.AlignCenter)
        horizontal_layout_main.addWidget(self.image_label)
        horizontal_layout_main.setStretch(0, 1)
        horizontal_layout_main.setStretch(1, 2)
        # =======================================================================================================
        # ---------panel s tlačítky-----------
        horizontal_layout_buttons = QHBoxLayout()
        # --------------Přidání tlačítek----------------
        self.start_button = QPushButton("Start")  # Vytvoření tlačítka start
        self.start_button.setStyleSheet("background-color: green;font-weight: bold; font-size: 12px;")
        horizontal_layout_buttons.addWidget(self.start_button)
        self.pause_button = QPushButton("Přeruš")  # Vytvoření tlačítka pause
        horizontal_layout_buttons.addWidget(self.pause_button)
        self.resume_button = QPushButton("Návrat")  # Vytvoření tlačítka resume
        horizontal_layout_buttons.addWidget(self.resume_button)
        self.stop_button = QPushButton("Stop")  # Vytvoření tlačítka stop
        self.stop_button.setStyleSheet("background-color: red")
        horizontal_layout_buttons.addWidget(self.stop_button)
        vertical_layout.addLayout(horizontal_layout_main)
        vertical_layout.addLayout(horizontal_layout_buttons)
        # TODO přidejte do GUI progresbar a provádějte jeho aktualizaci - viz metoda update_progress()
        self.setCentralWidget(central_widget)

    # obslužné metody pro úpravu hodnot
    def update_radius(self):
        """
        Metoda pro aktualizaci poloměru kruhu
        :return:
        """
        self.R = self.radius_input.value()*0.02
        self.radius_label.setText(f"Poloměr: {self.R:.2f}")
        self.setting_updated=True
    def update_T(self):
        """
        Metoda pro aktualizaci doby simulace
        :return:
        """
        self.T = self.T_input.value()
        self.setting_updated=True
    def update_V0(self):
        """
        Metoda pro aktualizaci počáteční rychlosti
        :return:
        """
        self.V0 = self.V0_input.value()
        self.setting_updated=True
    def update_nu(self):
        """
        Metoda pro aktualizaci viskozity
        :return:
        """
        self.nu = self.nu_input.value()
        self.setting_updated=True
    # Sloty pro ovládání výpočtu
    def start_computation(self):
        """
        Spustí výpočet
        :return:
        """
        if self.setting_updated: # pokud byly změny v GUI, tak je potřeba vytvořit nový worker
            self.worker = ComputationalLBM(T=self.T, u0=self.V0, nu=self.nu, R=self.R)  # Vytvoření instance výpočetního vlákna
            self.worker.progress_updated.connect(self.update_progress)  # Připojení slotu k signálu - update prograsu
            self.worker.velocity_field_updated.connect(self.update_image)  # Připojení slotu k signálu - update hodnoty
            self.worker.finished.connect(self.computation_finished)  # Připojení slotu k signálu - konec výpočtu

        self.statusBar().showMessage("Výpočet probíhá...") # Zobrazíme zprávu
        self.worker.current_time = 0 # vynuluj čas workeru
        if not self.worker.isRunning():  # Spustíme vlákno pouze, pokud neběží
            self.worker.start()
        self.setting_updated = False  # nastavíme, že změny byly zpracovány
        # TODO Bylo by dobré v průběhu výpočtu zakázat změny parametrů výpočtu - udělejte to
        # TODO Přidejte logiku pro zakázání-povolení stisknutí některých tlačítek - když není výpočet spuštěn, tak je zakážete

    def pause_computation(self):
        """
        Pozastaví výpočet
        :return:
        """
        self.worker.pause() # volej metodu pause
        self.statusBar().showMessage("Výpočet pozastaven...")

    def resume_computation(self):
        """
        Obnoví výpočet
        :return:
        """
        self.worker.resume()
        self.statusBar().showMessage("Výpočet obnoven...")

    def stop_computation(self):
        """
        Přeruší výpočet
        :return:
        """
        self.worker.stop()
        self.statusBar().showMessage("Výpočet přerušen...")

    @pyqtSlot(float)
    def update_progress(self, act_time):
        """
        Metoda pro aktualizaci postupu výpočtu
        :param act_time:
        :return:
        """
        progress = int((act_time / self.T) * 100)
        self.statusBar().showMessage(f"Progres: {progress}% - čas: {act_time:.2f}")
        # TODO:přidejte do GUI progresbar a provádějte jeho aktualizaci

    @pyqtSlot(QImage)
    def update_image(self, image):
        """
        Metoda pro aktualizaci obrázku
        :param image: vstupní vypočtený obrázek QImage
        :return:
        """
        pixmap = QPixmap.fromImage(image) # převedení QImage na QPixmap
        scaled_pixmap = pixmap.scaled(self.image_label.size(), Qt.KeepAspectRatio) # zvětšení obrázku dle velikosti labelu
        self.image_label.setPixmap(scaled_pixmap) # nastavení obrázku do labelu

    @pyqtSlot()
    def computation_finished(self):
        """
        Obsluha toho, co se má stát po dokončení výpočtu
        :return:
        """
        self.statusBar().showMessage("Výpočet dokončen!")
        # TODO povolte zakázaná tlačítka



if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = SimulationGui()
    window.show()
    sys.exit(app.exec_())
