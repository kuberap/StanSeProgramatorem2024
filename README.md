# Staň se na den programátorem


Tento projekt obsahuje různé výpočetní simulace implementované v Pythonu pomocí knihovny PyQt5 pro tvorbu GUI a dalších podpůrných knihoven. 
Cílem je ukázat, jak můžete využít Python k vytvoření jednoduchých výpočetních simulací a vizualizací pomocí GUI.

## Obsah

- **Monte Carlo Simulace pro π:** 
  - Jednoduchá GUI aplikace, která odhaduje hodnotu čísla π pomocí Monte Carlo metody.
  - Ukázka použití vláken pro výpočetní úlohy.
  
- **Proudění pomocí LBM:**
  - Komplexnější úloha pro simulaci proudění, který používá LBM. Obsahuje ovládací prvky pro úpravu parametrů simulace a zobrazení výsledků.
  - Používáme objektový princip zapouzdření a oddělujeme výpočet od GUI.

## Instalace
0. Vstupte do složky, kde chcete projekt uložit.
1. Stáhněte si tento repozitář pomocí příkazu:`git clone https://github.com/kuberap/StanSeProgramatorem2024`
a vstoupíme do složky: `cd StanSeProgramatorem2024/`
2. Vytvořte virtuální prostředí pomocí: `python -m venv venv`
3. Aktivujte virtuální prostředí pomocí: `venv\Scripts\activate` (Windows) nebo `source venv/bin/activate` (Linux)
a nainstalujte potřebné knihovny pomocí: `pip install -r requirements.txt`
## Spuštění
1. Změňte adresář `cd Src`
2. V adresáři jsou následující zdrojové soubory:

| **Soubor**               | **Popis**                                         | **Info**         |
|--------------------------|--------------------------------------------------|--------------------|
| `aplikace.py`            | Hlavní soubor aplikace, který obsahuje GUI.      | komplexní - *RUN*    |
| `computationalLBM.py`    | Adaptér pro LBM s použitím vláken   | pomocná            |
| `LBM.py`                 | Implementace LBM (Lattice Boltzmann Method).     | pomocná          |
| `vypocetPI.py`           | Výpočet čísla π pomocí Monte Carlo.              | snadná - *RUN*           |
| `vypocetPI_vlakna.py`    | Paralelní výpočet čísla π s využitím vláken.     | snadná - *RUN*           |

3. Spusťte aplikaci pomocí: `python vypocetPI.py` (nebo jiná aplikace s příznakem *RUN*)

