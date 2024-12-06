import sys
import random
from PyQt5.QtWidgets import QApplication, QLabel, QPushButton, QVBoxLayout, QWidget


class MainWindow(QWidget):
    """Hlavní GUI okno."""
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Výpočet \u03c0 ") # Nastavení titulku okna
        self.resize(300, 150)

        # Vytvoření tlačítka a labelu
        self.label = QLabel("Stiskněte 'Start' pro výpočet \u03c0", self)
        self.label.setStyleSheet("font-size: 14px;")
        self.start_button = QPushButton("Start", self)
        self.start_button.setStyleSheet("background-color: green;font-weight: bold; font-size: 16px;")
        self.start_button.clicked.connect(self.calculate_pi)

        # Rozvržení
        layout = QVBoxLayout(self)
        layout.addWidget(self.label)
        layout.addWidget(self.start_button)

    def calculate_pi(self):
        """Výpočet čísla pi přímo v hlavním vlákně."""
        self.label.setText("Výpočet probíhá...")
        self.start_button.setEnabled(False)  # Deaktivuje tlačítko během výpočtu

        inside = 0
        total_points = 10000000

        for i in range(1, total_points + 1):
            x = random.random()
            y = random.random()
            if x**2 + y**2 <= 1:
                inside += 1

            # aktualizace každých deset procent iterecí
            if i % (total_points // 10) == 0:
                pi_estimate = 4 * inside / i # spočti odhad čísla PI
                self.label.setText(f"Iterace: {i}, \u03c0 ≈ {pi_estimate:.6f}") # aktualizuj label
                # Pokud tu nebude následující řádek, GUI se zasekne
                # TODO vyzkoušejte zakomentovat a sledujte chování GUI
                QApplication.processEvents()  # Občasné zpracování GUI událostí

        # generujeme konečný odhad
        pi_final = 4 * inside / total_points
        self.label.setText(f"Odhad \u03c0 ≈ {pi_final:.6f}")
        self.start_button.setEnabled(True)


if __name__ == "__main__":
    # Vytvoří a spustí hlavní GUI aplikaci
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
